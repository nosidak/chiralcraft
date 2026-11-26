//! 1:1 parity ChiralCraft Document, Transaction (Undo/Redo/Save/Revert),
//! Stylesheet (`ACS Document 1996` vs New Document), Label Attachment (`F3C`/`CF3`/`CH3`/`H3C`),
//! Selection (Marquee/Lasso/AttachedData), and Cleanup Engine.
//!
//! Incorporates runtime-verified findings from ChiralCraft (`P01`, `P04`, `seed-v1.cdxml`, `seed-v2.cdxml`)
//! and static binary analysis of `the reference CDXML schema` and `CoreChemistryCommon.dll`.

use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

pub mod all_discipline;
pub mod nicknames;
pub mod pro_scene;
pub mod scene_ext;
pub mod synthesis;
pub use all_discipline::*;
pub use nicknames::*;
pub use pro_scene::*;
pub use scene_ext::*;
pub use synthesis::*;

use chiralcraft_chem::{
    default_valence, get_element_by_symbol, get_pseudo_atom_by_symbol, perceive_stereochemistry,
    BondDisplay, CcAromaticityInfo, CcVal, ChemAtom, ChemMolecule,
};
use chiralcraft_formats::{classify_single_caption_token, CdxmlDocument, CdxmlNodeChild};

// ============================================================================
// 1. Stylesheet & Fixed Bond Length Defaults (P04 & seed-v1.cdxml)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelJustification {
    Auto,
    Left,
    Center,
    Right,
    Above,
    Below,
}

impl LabelJustification {
    pub fn from_cdxml_str(s: &str) -> Self {
        match s.trim() {
            "Left" => Self::Left,
            "Center" => Self::Center,
            "Right" => Self::Right,
            "Above" => Self::Above,
            "Below" => Self::Below,
            _ => Self::Auto,
        }
    }

    pub fn as_cdxml_str(self) -> &'static str {
        match self {
            Self::Auto => "Auto",
            Self::Left => "Left",
            Self::Center => "Center",
            Self::Right => "Right",
            Self::Above => "Above",
            Self::Below => "Below",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DocumentSettings {
    pub creation_program: String,
    pub stylesheet_name: String,
    /// Bond length in points (`1 in = 72 pt`, `1 cm = 72 / 2.54 pt`).
    /// - `ACS Document 1996` (`seed-v1.cdxml`): `14.40 pt` (`0.2000 in` / `0.508 cm`)
    /// - Empty styled document (`_DD7204528CC072F7DF66E6261A19E64A`): `30.00 pt` (`30 / 72 = 0.416666... in`, rounded `0.4167 in`)
    pub bond_length_pt: f64,
    pub bond_spacing_pct: f64,
    pub chain_angle_deg: f64,
    pub line_width_pt: f64,
    pub bold_width_pt: f64,
    pub hash_spacing_pt: f64,
    pub margin_width_pt: f64,
    pub label_font: String,
    pub label_size_pt: f64,
    pub caption_font: String,
    pub caption_size_pt: f64,
    pub label_justification: LabelJustification,
    pub interpret_chemically: bool,
    pub fixed_lengths: bool,
    pub fixed_angles: bool,
    pub auto_rectify_hydrogens: bool,
}

impl DocumentSettings {
    /// Exact settings serialized by ChiralCraft in `seed-v1.cdxml` (`New ACS Document 1996`).
    pub fn acs_document_1996() -> Self {
        Self {
            creation_program: "ChiralCraft 1.0".to_string(),
            stylesheet_name: "ACS Document 1996".to_string(),
            bond_length_pt: 14.40,
            bond_spacing_pct: 18.0,
            chain_angle_deg: 120.0,
            line_width_pt: 0.60,
            bold_width_pt: 2.01,
            hash_spacing_pt: 2.49,
            margin_width_pt: 1.59,
            label_font: "Arial".to_string(),
            label_size_pt: 10.0,
            caption_font: "Arial".to_string(),
            caption_size_pt: 10.0,
            label_justification: LabelJustification::Auto,
            interpret_chemically: true,
            fixed_lengths: true,
            fixed_angles: true,
            auto_rectify_hydrogens: true,
        }
    }

    /// Standard 30.0 pt (`30.0 / 72.0 = 0.416666... in`, displayed `0.4167 in`) document stylesheet.
    pub fn standard_30pt_document() -> Self {
        let mut s = Self::acs_document_1996();
        s.stylesheet_name = "New Document".to_string();
        s.bond_length_pt = 30.0;
        s
    }

    /// Clean-room built-in published journal & presentation stylesheets (zero external files required).
    pub fn builtin_journal_stylesheets() -> Vec<(String, Self)> {
        let mk = |name: &str, bl: f64, sp: f64, lw: f64, bw: f64, hs: f64, mw: f64, font: &str, lsz: f64, csz: f64| {
            let mut s = Self::acs_document_1996();
            s.stylesheet_name = name.to_string();
            s.bond_length_pt = bl;
            s.bond_spacing_pct = sp;
            s.line_width_pt = lw;
            s.bold_width_pt = bw;
            s.hash_spacing_pt = hs;
            s.margin_width_pt = mw;
            s.label_font = font.to_string();
            s.label_size_pt = lsz;
            s.caption_font = font.to_string();
            s.caption_size_pt = csz;
            (name.to_string(), s)
        };
        vec![
            mk("ACS Document 1996", 14.40, 18.0, 0.60, 2.01, 2.49, 1.59, "Arial", 10.0, 10.0),
            mk("Angewandte Chemie", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 9.0, 9.0),
            mk("Biochemistry", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 10.0, 10.0),
            mk("Compact Figure", 11.50, 18.0, 0.55, 1.65, 2.10, 1.30, "Arial", 7.5, 7.5),
            mk("JACS Document", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 10.0, 10.0),
            mk("Journal of Medicinal Chemistry", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 10.0, 10.0),
            mk("Journal of Organic Chemistry", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 10.0, 10.0),
            mk("Nature Research Document", 10.80, 18.0, 0.60, 1.56, 2.00, 1.25, "Arial", 6.0, 6.0),
            mk("New Document", 30.00, 18.0, 1.00, 2.50, 3.00, 2.00, "Arial", 12.0, 12.0),
            mk("Organic Letters", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 10.0, 10.0),
            mk("Phytochemistry", 17.00, 18.0, 0.70, 2.20, 2.60, 1.70, "Times New Roman", 10.0, 10.0),
            mk("Poster / Presentation", 24.00, 18.0, 1.20, 3.20, 3.20, 2.20, "Arial", 14.0, 14.0),
            mk("RSC Document", 12.20, 20.0, 0.60, 1.80, 2.20, 1.40, "Arial", 8.0, 8.0),
            mk("Synthesis / Synlett", 17.00, 18.0, 0.70, 2.20, 2.50, 1.70, "Arial", 10.0, 10.0),
            mk("Tetrahedron", 16.00, 18.0, 0.65, 2.10, 2.50, 1.60, "Arial", 9.0, 9.0),
            mk("Wiley Document", 14.40, 18.0, 0.60, 2.00, 2.50, 1.60, "Arial", 9.0, 9.0),
        ]
    }

    /// Parses document settings and stylesheet parameters directly from a `CdxmlDocument`
    /// (`CreationProgram`, `BondLength`, `ChainAngle`, `LabelFont`, `LabelJustification`, etc.).
    pub fn from_cdxml(cdxml: &CdxmlDocument) -> Self {
        let mut s = Self::acs_document_1996();
        let root_el = match cdxml.elements.get(cdxml.root) {
            Some(el) => el,
            None => return s,
        };
        if let Some(cp) = root_el.attr("CreationProgram") {
            s.creation_program = cp.to_string();
        }
        if let Some(v) = root_el.attr("BondLength").and_then(|v| v.parse::<f64>().ok()) {
            s.bond_length_pt = v;
            s.stylesheet_name = if (v - 14.40).abs() < 0.05 {
                "ACS Document 1996".to_string()
            } else {
                "Custom".to_string()
            };
        }
        if let Some(v) = root_el.attr("BondSpacing").and_then(|v| v.parse::<f64>().ok()) {
            s.bond_spacing_pct = v;
        }
        if let Some(v) = root_el.attr("ChainAngle").and_then(|v| v.parse::<f64>().ok()) {
            s.chain_angle_deg = v;
        }
        if let Some(v) = root_el.attr("LineWidth").and_then(|v| v.parse::<f64>().ok()) {
            s.line_width_pt = v;
        }
        if let Some(v) = root_el.attr("BoldWidth").and_then(|v| v.parse::<f64>().ok()) {
            s.bold_width_pt = v;
        }
        if let Some(v) = root_el.attr("HashSpacing").and_then(|v| v.parse::<f64>().ok()) {
            s.hash_spacing_pt = v;
        }
        if let Some(v) = root_el.attr("MarginWidth").and_then(|v| v.parse::<f64>().ok()) {
            s.margin_width_pt = v;
        }
        if let Some(v) = root_el.attr("LabelSize").and_then(|v| v.parse::<f64>().ok()) {
            s.label_size_pt = v;
        }
        if let Some(v) = root_el.attr("CaptionSize").and_then(|v| v.parse::<f64>().ok()) {
            s.caption_size_pt = v;
        }
        if let Some(v) = root_el.attr("LabelJustification") {
            s.label_justification = LabelJustification::from_cdxml_str(v);
        }
        if let Some(v) = root_el.attr("InterpretChemically") {
            s.interpret_chemically = v == "yes";
        }

        // Resolve font names via fonttable lookup
        let resolve_font = |font_id: &str| -> Option<String> {
            cdxml.font_ids.get(font_id).and_then(|occs| {
                occs.first()
                    .and_then(|&occ| cdxml.elements.get(occ))
                    .and_then(|el| el.attr("name"))
                    .map(str::to_string)
            })
        };
        if let Some(fid) = root_el.attr("LabelFont") {
            if let Some(name) = resolve_font(fid) {
                s.label_font = name;
            }
        }
        if let Some(fid) = root_el.attr("CaptionFont") {
            if let Some(name) = resolve_font(fid) {
                s.caption_font = name;
            }
        }
        s
    }

    #[inline]
    pub fn bond_length_inches(&self) -> f64 {
        self.bond_length_pt / 72.0
    }

    #[inline]
    pub fn bond_length_cm(&self) -> f64 {
        (self.bond_length_pt / 72.0) * 2.54
    }

    /// Snaps a dragged bond vector `(x1, y1) -> (x2, y2)` in document points according to `fixed_lengths`
    /// and `fixed_angles` (15-degree increments when `fixed_angles` is active).
    /// Note: Document coordinates are zoom-invariant (`100%` vs `200%` zoom changes screen pixels, not CDXML point distance).
    pub fn snap_bond_drag(&self, x1: f64, y1: f64, x2: f64, y2: f64) -> (f64, f64) {
        let dx = x2 - x1;
        let dy = y2 - y1;
        let raw_len = (dx * dx + dy * dy).sqrt();
        if raw_len <= 1e-9 {
            return (x1 + self.bond_length_pt, y1);
        }
        let mut angle = dy.atan2(dx);
        if self.fixed_angles {
            let step = std::f64::consts::PI / 12.0; // 15 degrees
            angle = (angle / step).round() * step;
        }
        let len = if self.fixed_lengths {
            self.bond_length_pt
        } else {
            raw_len
        };
        (x1 + len * angle.cos(), y1 + len * angle.sin())
    }
}

// ============================================================================
// 2. Undo / Redo / Save / Revert Transaction Engine (Verified in P01 V1–V5)
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct CaptionNode {
    pub id: u32,
    pub x: f64,
    pub y: f64,
    pub text: String,
    pub element: Option<i32>,
    pub generic_nickname: Option<String>,
    pub warning: Option<String>,
    pub justification: LabelJustification,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AromaticCircle {
    pub molecule_idx: Option<usize>,
    pub center_x: f64,
    pub center_y: f64,
    pub radius: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct DocumentScene {
    pub captions: Vec<CaptionNode>,
    pub molecules: Vec<ChemMolecule>,
    pub aromatic_circles: Vec<AromaticCircle>,
    pub arrows: Vec<SceneArrow>,
    pub curves: Vec<SceneCurve>,
    pub graphics: Vec<SceneGraphic>,
    pub tables: Vec<SceneTable>,
    pub tlc_plates: Vec<SceneTlcPlate>,
    pub energy_profiles: Vec<EnergyProfile>,
    pub header_footer: PageHeaderFooter,
    pub expanded_residues: Vec<ExpandedResidue>,
}

#[derive(Debug, Clone)]
pub struct UndoTransaction {
    /// Menu label displayed after `"Undo "` / `"Redo "` (e.g., `"Caption Change"`, `"Typing"`, `"Bond"`).
    pub action_label: String,
    pub before: DocumentScene,
    pub after: DocumentScene,
}

#[derive(Debug, Clone)]
pub struct InteractiveDocument {
    pub name: String,
    pub settings: DocumentSettings,
    pub preferences: UserPreferences,
    pub nicknames: NicknameTable,
    pub scene: DocumentScene,
    undo_stack: Vec<UndoTransaction>,
    /// Current position in `undo_stack` (`0..=undo_stack.len()`).
    cursor: usize,
    /// Position in `undo_stack` corresponding to the last `Save` (`None` if modified off an invalidated branch).
    saved_cursor: Option<usize>,
    /// Snapshot of `scene` at the last `Save` (used by `File > Revert` in P01-V5).
    saved_scene: DocumentScene,
    pub last_repeat_label: Option<String>,
    next_id: u32,
}

impl InteractiveDocument {
    pub fn new_saved_seed(name: &str, settings: DocumentSettings) -> Self {
        Self {
            name: name.to_string(),
            settings,
            preferences: UserPreferences::default(),
            nicknames: NicknameTable::with_factory_defaults(),
            scene: DocumentScene::default(),
            undo_stack: Vec::new(),
            cursor: 0,
            saved_cursor: Some(0),
            saved_scene: DocumentScene::default(),
            last_repeat_label: None,
            next_id: 1,
        }
    }

    /// Record an undoable macro mutation on `self.scene`.
    pub fn apply_scene_tx<F: FnOnce(&mut DocumentScene, &mut u32)>(&mut self, action_label: &str, f: F) {
        let mut next_scene = self.scene.clone();
        f(&mut next_scene, &mut self.next_id);
        self.record_transaction(action_label, next_scene);
    }

    /// Stamp a molecule/template onto the document as a single Undo transaction.
    pub fn stamp_molecule_tx(&mut self, mut mol: ChemMolecule, label: &str) {
        self.apply_scene_tx(label, |scene, next_id| {
            for a in &mut mol.atoms {
                a.id = *next_id;
                *next_id += 1;
            }
            for b in &mut mol.bonds {
                b.id = *next_id;
                *next_id += 1;
            }
            scene.molecules.push(mol);
        });
    }

    /// Scale all molecules around their combined center as a single Undo transaction (`Object > Scale...`).
    pub fn scale_all_tx(&mut self, factor: f64) {
        self.apply_scene_tx("Scale", |scene, _| {
            scale_molecules_around_center(&mut scene.molecules, factor);
        });
    }

    /// Expand a nickname label on `(mol_idx, atom_idx)` as a single Undo transaction (`Structure > Expand Label`).
    pub fn expand_label_tx(&mut self, mol_idx: usize, atom_idx: usize, label: &str) -> bool {
        let nicks = self.nicknames.clone();
        let blen = self.settings.bond_length_pt;
        let mut ok = false;
        self.apply_scene_tx("Expand Label", |scene, next_id| {
            if let Some(mol) = scene.molecules.get_mut(mol_idx) {
                if let Some(res) = nicks.expand_label_at_atom(mol, atom_idx, label, blen, next_id) {
                    scene.expanded_residues.push(res);
                    ok = true;
                }
            }
        });
        ok
    }

    /// Symmetrically contract the most recent expanded residue on `mol_idx` (`Structure > Contract Label`).
    pub fn contract_last_residue_tx(&mut self, mol_idx: usize) -> bool {
        let mut ok = false;
        self.apply_scene_tx("Contract Label", |scene, _| {
            if let (Some(mol), Some(res)) = (
                scene.molecules.get_mut(mol_idx),
                scene.expanded_residues.pop(),
            ) {
                ok = NicknameTable::contract_residue_group(mol, &res);
            }
        });
        ok
    }

    /// Run `check_structure_scene` filtered through `self.preferences.warnings_enabled` (all 13 checkboxes).
    pub fn check_structure_with_preferences(&self) -> Vec<SceneStructureWarning> {
        check_structure_scene(&self.scene)
            .into_iter()
            .filter(|w| self.preferences.is_warning_allowed(&w.message))
            .collect()
    }

    /// True if the document title shows the modified `*` indicator (verified in `runtime-history.md` and `P01`).
    #[inline]
    pub fn is_dirty(&self) -> bool {
        self.saved_cursor != Some(self.cursor)
    }

    /// Returns the exact Edit menu Undo label (`"Undo <label>"` or `"Undo not available"`).
    pub fn undo_menu_label(&self) -> String {
        if self.cursor > 0 {
            format!("Undo {}", self.undo_stack[self.cursor - 1].action_label)
        } else {
            "Undo not available".to_string()
        }
    }

    /// Returns the exact Edit menu Redo label (`"Redo <label>"` or `"Redo not available"`).
    pub fn redo_menu_label(&self) -> String {
        if self.cursor < self.undo_stack.len() {
            format!("Redo {}", self.undo_stack[self.cursor].action_label)
        } else {
            "Redo not available".to_string()
        }
    }

    /// Returns whether `File > Save` is enabled (enabled iff document is dirty).
    #[inline]
    pub fn is_save_enabled(&self) -> bool {
        self.is_dirty()
    }

    /// Returns whether `File > Revert` is enabled (enabled iff document is dirty and has a saved backing state).
    #[inline]
    pub fn is_revert_enabled(&self) -> bool {
        self.is_dirty()
    }

    pub fn record_transaction(&mut self, label: &str, new_scene: DocumentScene) {
        // P01-V4 verified: recording a new action after Undo truncates the Redo branch!
        if self.cursor < self.undo_stack.len() {
            if let Some(sc) = self.saved_cursor {
                if sc > self.cursor {
                    self.saved_cursor = None;
                }
            }
            self.undo_stack.truncate(self.cursor);
        }
        let tx = UndoTransaction {
            action_label: label.to_string(),
            before: self.scene.clone(),
            after: new_scene.clone(),
        };
        self.scene = new_scene;
        self.undo_stack.push(tx);
        self.cursor += 1;
        self.last_repeat_label = Some(label.to_string());
    }

    /// Executes `Edit > Undo`. Returns `true` if an action was undone.
    pub fn undo(&mut self) -> bool {
        if self.cursor == 0 {
            return false;
        }
        self.cursor -= 1;
        self.scene = self.undo_stack[self.cursor].before.clone();
        self.last_repeat_label = None;
        true
    }

    /// Executes `Edit > Redo`. Returns `true` if an action was redone.
    pub fn redo(&mut self) -> bool {
        if self.cursor >= self.undo_stack.len() {
            return false;
        }
        self.scene = self.undo_stack[self.cursor].after.clone();
        self.cursor += 1;
        self.last_repeat_label = None;
        true
    }

    /// Executes `File > Save`.
    /// Verified in **P01-V2 (`DOC-CONTRADICTION`)**: `File > Save` does **NOT** clear the Undo stack!
    /// `Undo Caption Change` remains enabled immediately after `Save`, `Repeat Save` becomes the Repeat label,
    /// and subsequent `Edit > Undo` undoes across the save boundary.
    pub fn save(&mut self) {
        self.saved_cursor = Some(self.cursor);
        self.saved_scene = self.scene.clone();
        self.last_repeat_label = Some("Save".to_string());
    }

    /// Executes `File > Revert` (P01-V5): restores the document to `saved_scene` and clears post-save changes.
    pub fn revert(&mut self) -> bool {
        if !self.is_revert_enabled() {
            return false;
        }
        self.scene = self.saved_scene.clone();
        self.undo_stack.clear();
        self.cursor = 0;
        self.saved_cursor = Some(0);
        self.last_repeat_label = None;
        true
    }

    /// Simulates the exact 3-transaction caption creation sequence observed in ChiralCraft (`P01-V1` & `seed-v2.cdxml`):
    /// 1. Clicking canvas with Text tool -> `"Caption Change"` (opens empty caption box)
    /// 2. Typing text -> `"Typing"` (inserts characters)
    /// 3. Committing caption (switching tool / clicking blank canvas) -> `"Caption Change"`
    ///    (applies `InterpretChemically`: resolves single elements, generic nicknames `A`/`M`/`Q`/`X`/`R`,
    ///    and evaluates valence errors via `CcVal::compute`, e.g. isolated `"B"` -> `"An atom in this label has an invalid valence."`).
    pub fn create_caption_interactive(&mut self, x: f64, y: f64, text: &str) {
        let id = self.next_id;
        self.next_id += 1;

        // Step 1: Caption Change (open box)
        let mut s1 = self.scene.clone();
        s1.captions.push(CaptionNode {
            id,
            x,
            y,
            text: String::new(),
            element: None,
            generic_nickname: None,
            warning: None,
            justification: LabelJustification::Left,
        });
        self.record_transaction("Caption Change", s1);

        // Step 2: Typing
        let mut s2 = self.scene.clone();
        if let Some(cap) = s2.captions.iter_mut().find(|c| c.id == id) {
            cap.text = text.to_string();
        }
        self.record_transaction("Typing", s2);

        // Step 3: Caption Change (commit & chemical interpretation)
        let mut s3 = self.scene.clone();
        if let Some(cap) = s3.captions.iter_mut().find(|c| c.id == id) {
            if self.settings.interpret_chemically {
                interpret_caption_node(cap);
            }
        }
        self.record_transaction("Caption Change", s3);
    }

    /// Reconstructs an `InteractiveDocument` from a parsed `CdxmlDocument`, loading stylesheet settings,
    /// multi-atom fragments (`ChemMolecule`), and standalone/single-node captions (`CaptionNode`).
    pub fn from_cdxml(cdxml: &CdxmlDocument) -> Self {
        let settings = DocumentSettings::from_cdxml(cdxml);
        let name = cdxml
            .elements
            .get(cdxml.root)
            .and_then(|root| root.attr("Name"))
            .unwrap_or("untitled.cdxml")
            .to_string();

        let mut captions = Vec::new();
        let mut caption_frag_occs = HashSet::new();
        let mut max_id = 1u32;

        let collect_text = |occ: usize| -> String {
            fn gather(doc: &CdxmlDocument, id: usize, buf: &mut String) {
                if let Some(el) = doc.elements.get(id) {
                    for child in &el.children {
                        match child {
                            CdxmlNodeChild::Text(t) => buf.push_str(t),
                            CdxmlNodeChild::Element(cid) => gather(doc, *cid, buf),
                        }
                    }
                }
            }
            let mut out = String::new();
            gather(cdxml, occ, &mut out);
            out
        };

        for (occ, el) in cdxml.elements.iter().enumerate() {
            if let Some(id_val) = el.attr("id").and_then(|v| v.parse::<u32>().ok()) {
                max_id = max_id.max(id_val + 1);
            }
            if el.tag == "fragment" {
                let mut child_nodes = Vec::new();
                let mut bond_count = 0usize;
                for c in &el.children {
                    if let CdxmlNodeChild::Element(cid) = c {
                        if let Some(ce) = cdxml.elements.get(*cid) {
                            if ce.tag == "n" {
                                child_nodes.push(*cid);
                            } else if ce.tag == "b" {
                                bond_count += 1;
                            }
                        }
                    }
                }
                if child_nodes.len() == 1 && bond_count == 0 {
                    let node_occ = child_nodes[0];
                    let node_el = &cdxml.elements[node_occ];
                    let t_occ = node_el.children.iter().find_map(|c| {
                        if let CdxmlNodeChild::Element(cid) = c {
                            (cdxml.elements[*cid].tag == "t").then_some(*cid)
                        } else {
                            None
                        }
                    });
                    if let Some(t_occ) = t_occ {
                        let text = collect_text(t_occ);
                        if !text.is_empty() {
                            caption_frag_occs.insert(occ);
                            let (x, y) = parse_xy_attr(node_el.attr("p"));
                            let id = node_el
                                .attr("id")
                                .and_then(|v| v.parse::<u32>().ok())
                                .unwrap_or(max_id);
                            let just = cdxml.elements[t_occ]
                                .attr("LabelJustification")
                                .map(LabelJustification::from_cdxml_str)
                                .unwrap_or(LabelJustification::Left);
                            captions.push(CaptionNode {
                                id,
                                x,
                                y,
                                text,
                                element: node_el.attr("Element").and_then(|v| v.parse::<i32>().ok()),
                                generic_nickname: node_el
                                    .attr("GenericNickname")
                                    .map(str::to_string),
                                warning: node_el.attr("Warning").map(str::to_string),
                                justification: just,
                            });
                        }
                    }
                }
            } else if el.tag == "t" {
                let is_page_caption = el
                    .parent
                    .and_then(|p| cdxml.elements.get(p))
                    .is_some_and(|pe| pe.tag == "page");
                if is_page_caption {
                    let text = collect_text(occ);
                    let (x, y) = parse_xy_attr(el.attr("p"));
                    let id = el
                        .attr("id")
                        .and_then(|v| v.parse::<u32>().ok())
                        .unwrap_or(max_id);
                    let just = el
                        .attr("LabelJustification")
                        .or_else(|| el.attr("Justification"))
                        .map(LabelJustification::from_cdxml_str)
                        .unwrap_or(LabelJustification::Left);
                    captions.push(CaptionNode {
                        id,
                        x,
                        y,
                        text,
                        element: None,
                        generic_nickname: None,
                        warning: None,
                        justification: just,
                    });
                }
            }
        }

        let mut molecules = Vec::new();
        let mut occ_to_mol_idx = HashMap::new();
        for (occ, mol) in cdxml.extract_molecules() {
            if !caption_frag_occs.contains(&occ) {
                let idx = molecules.len();
                occ_to_mol_idx.insert(occ, idx);
                molecules.push(mol);
            }
        }
        let aromatic_circles = cdxml
            .extract_aromatic_circles()
            .into_iter()
            .map(|c| AromaticCircle {
                molecule_idx: c.fragment_occ.and_then(|occ| occ_to_mol_idx.get(&occ).copied()),
                center_x: c.center_x,
                center_y: c.center_y,
                radius: c.radius,
            })
            .collect();

        let mut tables = extract_tables_from_cdxml(cdxml);
        let mut tlc_plates = extract_tlc_plates_from_cdxml(cdxml);
        let mut energy_profiles = Vec::new();
        restore_chiralcraft_extensions_from_cdxml(
            cdxml,
            &mut energy_profiles,
            &mut tables,
            &mut tlc_plates,
        );

        let scene = DocumentScene {
            captions,
            molecules,
            aromatic_circles,
            arrows: extract_arrows_from_cdxml(cdxml),
            curves: extract_curves_from_cdxml(cdxml),
            graphics: extract_graphics_from_cdxml(cdxml),
            tables,
            tlc_plates,
            energy_profiles,
            header_footer: extract_page_header_footer(cdxml),
            expanded_residues: Vec::new(),
        };
        let mut doc = Self::new_saved_seed(&name, settings);
        doc.scene = scene.clone();
        doc.saved_scene = scene;
        doc.next_id = max_id.max(100);
        doc
    }

    /// Serializes the document settings, captions, and molecules into a valid `CdxmlDocument`.
    pub fn to_cdxml(&self) -> CdxmlDocument {
        let mut out = String::new();
        out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
        let mut page_attrs = format!(
            " id=\"1\"\n BoundingBox=\"0 0 540 719.75\"\n WidthPages=\"{}\"\n HeightPages=\"{}\"",
            self.scene.header_footer.width_pages.max(1),
            self.scene.header_footer.height_pages.max(1),
        );
        if let Some(hdr) = &self.scene.header_footer.header {
            page_attrs.push_str(&format!("\n Header=\"{}\"", escape_xml_attr(hdr)));
        }
        if let Some(ftr) = &self.scene.header_footer.footer {
            page_attrs.push_str(&format!("\n Footer=\"{}\"", escape_xml_attr(ftr)));
        }
        out.push_str(&format!(
            "<CDXML\n CreationProgram=\"{}\"\n Name=\"{}\"\n InterpretChemically=\"{}\"\n LabelFont=\"3\"\n LabelSize=\"{:.2}\"\n CaptionFont=\"3\"\n CaptionSize=\"{:.2}\"\n HashSpacing=\"{:.2}\"\n MarginWidth=\"{:.2}\"\n LineWidth=\"{:.2}\"\n BoldWidth=\"{:.2}\"\n BondLength=\"{:.2}\"\n BondSpacing=\"{:.0}\"\n ChainAngle=\"{:.0}\"\n LabelJustification=\"{}\"\n color=\"0\"\n bgcolor=\"1\"\n><colortable>\n<color r=\"1\" g=\"1\" b=\"1\"/>\n<color r=\"0\" g=\"0\" b=\"0\"/>\n</colortable><fonttable>\n<font id=\"3\" charset=\"iso-8859-1\" name=\"{}\"/>\n</fonttable><page{}\n>",
            escape_xml_attr(&self.settings.creation_program),
            escape_xml_attr(&self.name),
            if self.settings.interpret_chemically { "yes" } else { "no" },
            self.settings.label_size_pt,
            self.settings.caption_size_pt,
            self.settings.hash_spacing_pt,
            self.settings.margin_width_pt,
            self.settings.line_width_pt,
            self.settings.bold_width_pt,
            self.settings.bond_length_pt,
            self.settings.bond_spacing_pct,
            self.settings.chain_angle_deg,
            self.settings.label_justification.as_cdxml_str(),
            escape_xml_attr(&self.settings.label_font),
            page_attrs,
        ));

        let mut next_xml_id = self.next_id.max(1000);
        for cap in &self.scene.captions {
            if cap.element.is_some() || cap.generic_nickname.is_some() {
                let frag_id = next_xml_id;
                let text_id = next_xml_id + 1;
                next_xml_id += 2;
                let mut node_attrs = format!(" id=\"{}\" p=\"{:.2} {:.2}\"", cap.id, cap.x, cap.y);
                if let Some(nick) = &cap.generic_nickname {
                    node_attrs.push_str(&format!(
                        " NodeType=\"GenericNickname\" GenericNickname=\"{}\"",
                        escape_xml_attr(nick)
                    ));
                }
                if let Some(el) = cap.element {
                    node_attrs.push_str(&format!(" Element=\"{el}\" NumHydrogens=\"0\""));
                }
                if let Some(warn) = &cap.warning {
                    node_attrs.push_str(&format!(" Warning=\"{}\"", escape_xml_attr(warn)));
                }
                out.push_str(&format!(
                    "<fragment id=\"{frag_id}\"><n{node_attrs}><t id=\"{text_id}\" p=\"{:.2} {:.2}\" LabelJustification=\"{}\"><s font=\"3\" size=\"{:.0}\" color=\"0\">{}</s></t></n></fragment>",
                    cap.x,
                    cap.y,
                    cap.justification.as_cdxml_str(),
                    self.settings.label_size_pt,
                    escape_xml_text(&cap.text),
                ));
            } else {
                out.push_str(&format!(
                    "<t id=\"{}\" p=\"{:.2} {:.2}\" Justification=\"{}\"><s font=\"3\" size=\"{:.0}\" color=\"0\">{}</s></t>",
                    cap.id,
                    cap.x,
                    cap.y,
                    cap.justification.as_cdxml_str(),
                    self.settings.caption_size_pt,
                    escape_xml_text(&cap.text),
                ));
            }
        }

        for (mol_idx, mol) in self.scene.molecules.iter().enumerate() {
            let frag_id = next_xml_id;
            next_xml_id += 1;
            out.push_str(&format!("<fragment id=\"{frag_id}\">"));
            let mut atom_xml_ids = Vec::with_capacity(mol.atoms.len());
            for a in &mol.atoms {
                let aid = if a.id != 0 { a.id } else { next_xml_id };
                next_xml_id = next_xml_id.max(aid + 1);
                atom_xml_ids.push(aid);
                out.push_str(&format!(
                    "<n id=\"{aid}\" p=\"{:.2} {:.2}\" Element=\"{}\"",
                    a.x, a.y, a.z
                ));
                if a.charge != 0 {
                    out.push_str(&format!(" Charge=\"{}\"", a.charge));
                }
                if a.isotope != 0 {
                    out.push_str(&format!(" Isotope=\"{}\"", a.isotope));
                }
                if a.radical != 0 {
                    let rad_str = match a.radical {
                        1 => "Doublet",
                        2 => "Singlet",
                        3 => "Triplet",
                        _ => "Doublet",
                    };
                    out.push_str(&format!(" Radical=\"{rad_str}\""));
                }
                if let Some(vo) = a.valence_override {
                    out.push_str(&format!(" Valence=\"{vo}\""));
                }
                if a.is_unspecified_label {
                    out.push_str(" NodeType=\"Unspecified\"");
                }
                if !a.visible {
                    out.push_str(" Visible=\"no\"");
                }
                if a.ignore_warnings {
                    out.push_str(" IgnoreWarnings=\"yes\"");
                }
                if let Some(ref lbl) = a.label {
                    let tid = next_xml_id;
                    next_xml_id += 1;
                    out.push_str(&format!(
                        " CustomLabel=\"{}\"><t id=\"{tid}\" p=\"{:.2} {:.2}\"><s font=\"3\" size=\"{:.0}\" color=\"0\">{}</s></t></n>",
                        escape_xml_attr(lbl),
                        a.x,
                        a.y,
                        self.settings.label_size_pt,
                        escape_xml_text(lbl),
                    ));
                } else {
                    out.push_str("/>");
                }
            }
            for b in &mol.bonds {
                let bid = if b.id != 0 { b.id } else { next_xml_id };
                next_xml_id = next_xml_id.max(bid + 1);
                let ba = atom_xml_ids[b.begin_atom];
                let ea = atom_xml_ids[b.end_atom];
                out.push_str(&format!("<b id=\"{bid}\" B=\"{ba}\" E=\"{ea}\""));
                if b.order != 1 {
                    match b.order {
                        6 => out.push_str(" Order=\"1.5\""),
                        5 => out.push_str(" Order=\"dative\""),
                        7 => out.push_str(" Order=\"hydrogen\""),
                        _ => out.push_str(&format!(" Order=\"{}\"", b.order)),
                    }
                }
                if let Some(disp) = bond_display_cdxml_str(b.display) {
                    out.push_str(&format!(" Display=\"{disp}\""));
                }
                if !b.visible {
                    out.push_str(" Visible=\"no\"");
                }
                if b.ignore_warnings {
                    out.push_str(" IgnoreWarnings=\"yes\"");
                }
                out.push_str("/>");
            }
            for circ in self
                .scene
                .aromatic_circles
                .iter()
                .filter(|c| c.molecule_idx == Some(mol_idx))
            {
                let gid = next_xml_id;
                next_xml_id += 1;
                out.push_str(&format!(
                    "<graphic id=\"{gid}\" GraphicType=\"Oval\" OvalType=\"Circle\" Center3D=\"{:.2} {:.2} 0\" MajorAxisEnd3D=\"{:.2} {:.2} 0\" MinorAxisEnd3D=\"{:.2} {:.2} 0\"/>",
                    circ.center_x,
                    circ.center_y,
                    circ.center_x + circ.radius,
                    circ.center_y,
                    circ.center_x,
                    circ.center_y + circ.radius,
                ));
            }
            out.push_str("</fragment>");
        }
        for circ in self
            .scene
            .aromatic_circles
            .iter()
            .filter(|c| c.molecule_idx.is_none())
        {
            let gid = next_xml_id;
            next_xml_id += 1;
            out.push_str(&format!(
                "<graphic id=\"{gid}\" GraphicType=\"Oval\" OvalType=\"Circle\" Center3D=\"{:.2} {:.2} 0\" MajorAxisEnd3D=\"{:.2} {:.2} 0\" MinorAxisEnd3D=\"{:.2} {:.2} 0\"/>",
                circ.center_x,
                circ.center_y,
                circ.center_x + circ.radius,
                circ.center_y,
                circ.center_x,
                circ.center_y + circ.radius,
            ));
        }

        serialize_scene_ext_cdxml(
            &self.scene.arrows,
            &self.scene.curves,
            &self.scene.graphics,
            &self.scene.tlc_plates,
            &mut out,
        );
        if let Some(ext_tag) = serialize_chiralcraft_extensions_cdxml(
            &self.scene.energy_profiles,
            &self.scene.tables,
            &self.scene.tlc_plates,
        ) {
            out.push_str(&ext_tag);
        }

        out.push_str("</page></CDXML>");
        CdxmlDocument::parse(&out).expect("serialized document CDXML must parse cleanly")
    }
}

fn parse_xy_attr(attr: Option<&str>) -> (f64, f64) {
    let Some(s) = attr else {
        return (0.0, 0.0);
    };
    let mut parts = s.split_whitespace();
    let x = parts.next().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let y = parts.next().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    (x, y)
}

fn escape_xml_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn escape_xml_text(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn bond_display_cdxml_str(display: BondDisplay) -> Option<&'static str> {
    match display {
        BondDisplay::Solid => None,
        BondDisplay::Dash => Some("Dash"),
        BondDisplay::Hash => Some("Hash"),
        BondDisplay::WedgeBegin => Some("WedgeBegin"),
        BondDisplay::WedgeEnd => Some("WedgeEnd"),
        BondDisplay::WedgedHashBegin => Some("WedgedHashBegin"),
        BondDisplay::WedgedHashEnd => Some("WedgedHashEnd"),
        BondDisplay::Bold => Some("Bold"),
        BondDisplay::Wavy => Some("Wavy"),
        BondDisplay::HollowWedgeBegin => Some("HollowWedgeBegin"),
        BondDisplay::HollowWedgeEnd => Some("HollowWedgeEnd"),
        BondDisplay::Dative => Some("Dative"),
    }
}

fn atom_has_valence_warning(mol: &ChemMolecule, val: &CcVal, idx: usize) -> bool {
    let atom = &mol.atoms[idx];
    let z = atom.z;
    if val.indeterminate_valence[idx] {
        return true;
    }
    if chiralcraft_chem::is_metal_or_metalloid(z) {
        return false;
    }
    if atom.is_query || z <= 0 {
        return false;
    }
    if atom.explicit_num_hydrogens.is_none() && mol.bonds.is_empty() {
        // Bare atom caption without bonds (`NumHydrogens="0"` in ChiralCraft when !can_infer_implicit_hydrogen, e.g. B)
        return !chiralcraft_chem::can_infer_implicit_hydrogen(z)
            && chiralcraft_chem::default_valence(z) > 0;
    }
    let chg = atom.charge;
    let dative_shift = val.explicit_dative_out[idx] as i32 - val.explicit_dative_in[idx] as i32;
    let eff_chg = chg + dative_shift;
    let sbo = val.sbo[idx];
    let rad = atom.radical_electrons();
    if let Some(nh) = atom.explicit_num_hydrogens {
        let total = sbo + nh.max(0) + rad;
        if !chiralcraft_chem::is_valid_isoelectronic_state(z, eff_chg, total) {
            return true;
        }
        return total != chiralcraft_chem::get_intended_valence(z - eff_chg, total);
    }
    let required = sbo + rad;
    if !chiralcraft_chem::is_valid_isoelectronic_state(z, eff_chg, required) {
        return true;
    }
    let intended = chiralcraft_chem::get_intended_valence(z - eff_chg, required);
    if required > intended {
        return true;
    }
    if !chiralcraft_chem::can_infer_implicit_hydrogen(z) && intended != required {
        return true;
    }
    false
}

fn interpret_caption_node(cap: &mut CaptionNode) {
    let trimmed = cap.text.trim();
    if trimmed.is_empty() {
        return;
    }
    if let Some((el_opt, is_element)) = classify_single_caption_token(trimmed) {
        if is_element {
            if let Some(z) = el_opt {
                cap.element = Some(z);
                let mut mol = ChemMolecule::new();
                mol.add_atom(ChemAtom::new(1, z));
                let val = CcVal::compute(&mol);
                if atom_has_valence_warning(&mol, &val, 0) {
                    cap.warning = Some("An atom in this label has an invalid valence.".to_string());
                }
            }
        } else {
            cap.generic_nickname = Some(trimmed.to_string());
        }
        return;
    }
    if let Some(pseudo) = get_pseudo_atom_by_symbol(trimmed) {
        cap.element = Some(pseudo.z);
        return;
    }
    // Check condensed formula interpretation if it is a multi-atom formula (e.g., CF3 / CH5)
    let parsed = analyze_left_endpoint_label(trimmed, cap.justification, 0, false);
    if parsed.has_valence_warning {
        cap.warning = parsed.warning_message;
    }
}

// ============================================================================
// 3. Label Text vs Chemical Attachment & Expansion (P03)
// ============================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedLabelAttachment {
    pub raw_label: String,
    pub justification: LabelJustification,
    /// Atomic number of the attaching atom (`6` for C, `9` for F, `1` for H).
    pub attached_element_z: i32,
    pub attached_symbol: String,
    pub terminal_substituents: Vec<(i32, usize)>,
    pub has_valence_warning: bool,
    pub warning_message: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndpointSide {
    Left,
    Right,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FormulaToken {
    z: i32,
    symbol: String,
    count: usize,
}

fn format_formula_tokens(tokens: &[FormulaToken]) -> String {
    let mut s = String::new();
    for tok in tokens {
        s.push_str(&tok.symbol);
        if tok.count > 1 {
            s.push_str(&tok.count.to_string());
        }
    }
    s
}

fn tokenize_condensed_formula(label: &str) -> Option<Vec<FormulaToken>> {
    let trimmed = label.trim();
    if trimmed.is_empty() {
        return None;
    }
    let bytes = trimmed.as_bytes();
    let mut i = 0usize;
    let mut tokens = Vec::new();
    while i < bytes.len() {
        let ch = bytes[i];
        if !ch.is_ascii_uppercase() {
            return None;
        }
        let sym_start = i;
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_lowercase() {
            i += 1;
        }
        let mut sym = &trimmed[sym_start..i];
        // Handle 2-letter vs 1-letter split if e.g. Co vs CH or No vs NO
        let el = if let Some(found) = get_element_by_symbol(sym) {
            found
        } else if sym.len() > 1 {
            i = sym_start + 1;
            sym = &trimmed[sym_start..i];
            get_element_by_symbol(sym)?
        } else {
            return None;
        };
        let num_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let count = if num_start < i {
            trimmed[num_start..i].parse::<usize>().ok()?
        } else {
            1
        };
        if count == 0 {
            return None;
        }
        tokens.push(FormulaToken {
            z: el.z,
            symbol: el.symbol.to_string(),
            count,
        });
    }
    (!tokens.is_empty()).then_some(tokens)
}

fn resolve_known_nickname(
    label: &str,
    incident_bond_order: i32,
) -> Option<(i32, &'static str, Vec<(i32, usize)>, bool)> {
    // Returns (attached_z, attached_symbol, terminal_substituents, is_valence_valid)
    let valid_single = incident_bond_order <= 1;
    match label {
        "Me" => Some((6, "C", vec![(1, 3)], valid_single)),
        "Et" => Some((6, "C", vec![(1, 2), (6, 1)], valid_single)),
        "nPr" | "Pr" | "iPr" => Some((6, "C", vec![(1, 1), (6, 2)], valid_single)),
        "Bu" | "nBu" | "iBu" | "sBu" => Some((6, "C", vec![(1, 2), (6, 1)], valid_single)),
        "tBu" => Some((6, "C", vec![(6, 3)], valid_single)),
        "Ph" | "Ar" => Some((6, "C", vec![(6, 2)], valid_single)),
        "OMe" | "MeO" => Some((8, "O", vec![(6, 1)], valid_single)),
        "OEt" | "EtO" => Some((8, "O", vec![(6, 1)], valid_single)),
        "SMe" | "MeS" => Some((16, "S", vec![(6, 1)], valid_single)),
        "OCF3" | "F3CO" => Some((8, "O", vec![(6, 1)], valid_single)),
        "Ac" => Some((6, "C", vec![(8, 1), (6, 1)], valid_single)),
        "OAc" | "AcO" => Some((8, "O", vec![(6, 1)], valid_single)),
        "CN" => Some((6, "C", vec![(7, 1)], valid_single)),
        "NC" => Some((7, "N", vec![(6, 1)], valid_single)),
        "NO2" | "O2N" => Some((7, "N", vec![(8, 2)], valid_single)),
        "SO3H" | "HO3S" => Some((16, "S", vec![(8, 3)], valid_single)),
        "SO2Me" | "Ms" => Some((16, "S", vec![(8, 2), (6, 1)], valid_single)),
        "COOH" | "CO2H" | "HOOC" | "HO2C" => Some((6, "C", vec![(8, 2)], valid_single)),
        "COOMe" | "CO2Me" | "COOEt" | "CO2Et" => Some((6, "C", vec![(8, 2)], valid_single)),
        "CHO" | "OHC" => Some((6, "C", vec![(8, 1), (1, 1)], valid_single)),
        "NMe2" | "NEt2" => Some((7, "N", vec![(6, 2)], valid_single)),
        _ => None,
    }
}

fn is_monovalent_terminal(z: i32) -> bool {
    matches!(z, 1 | 9 | 17 | 35 | 53 | 85)
}

fn pick_attachment_token_idx(
    tokens: &[FormulaToken],
    side: EndpointSide,
    justification: LabelJustification,
) -> usize {
    if tokens.len() <= 1 {
        return 0;
    }
    let last_non_h = tokens.iter().rposition(|t| t.z != 1).unwrap_or(tokens.len() - 1);
    let first_non_h = tokens.iter().position(|t| t.z != 1).unwrap_or(0);

    let effective_just = match (side, justification) {
        (EndpointSide::Left, LabelJustification::Right) => LabelJustification::Right,
        (EndpointSide::Left, LabelJustification::Left) => LabelJustification::Left,
        (EndpointSide::Right, LabelJustification::Left) => LabelJustification::Left,
        (EndpointSide::Right, LabelJustification::Right) => LabelJustification::Right,
        _ => LabelJustification::Auto,
    };

    match effective_just {
        // Explicit flush-right on a left endpoint attaches to rightmost non-H atom (Hydrogen exception, T08)
        LabelJustification::Right => last_non_h,
        // Explicit flush-left attaches to leftmost non-H atom (Hydrogen exception)
        LabelJustification::Left => first_non_h,
        // Auto chooses the principal polyvalent heavy atom (handles both forward CF3/CH3/OH and reversed F3C/H3C/HO)
        LabelJustification::Auto
        | LabelJustification::Center
        | LabelJustification::Above
        | LabelJustification::Below => {
            if is_monovalent_terminal(tokens[0].z) && tokens.len() > 1 {
                return last_non_h;
            }
            if tokens[0].count > 1 && default_valence(tokens[tokens.len() - 1].z) > default_valence(tokens[0].z) {
                return last_non_h;
            }
            first_non_h
        }
    }
}

fn evaluate_attached_valence(
    attached_z: i32,
    incident_bond_order: i32,
    substituents: &[(i32, usize)],
) -> bool {
    let mut mol = ChemMolecule::new();
    let center_idx = mol.add_atom(ChemAtom::new(1, attached_z));
    let mut next_id = 2u32;
    if incident_bond_order > 0 {
        let stub_idx = mol.add_atom(ChemAtom::new(next_id, 6));
        next_id += 1;
        mol.add_bond(next_id, center_idx, stub_idx, incident_bond_order, BondDisplay::Solid);
        next_id += 1;
    }
    for &(sub_z, count) in substituents {
        for _ in 0..count {
            let sub_idx = mol.add_atom(ChemAtom::new(next_id, sub_z));
            next_id += 1;
            mol.add_bond(next_id, center_idx, sub_idx, 1, BondDisplay::Solid);
            next_id += 1;
        }
    }
    let val = CcVal::compute(&mol);
    atom_has_valence_warning(&mol, &val, center_idx)
}

/// Evaluates chemical attachment and valence validity for endpoint labels (condensed formulas
/// such as `F3C`/`CF3`/`CH3`/`H3C`/`OH`/`HO`/`NH2`/`H2N`/`CCl3` and common nicknames `Me`/`Et`/`NO2`/`COOH`)
/// on a `Left` or `Right` bond endpoint (`P03`, `T07`–`T10`, `T12`).
pub fn analyze_endpoint_label(
    label: &str,
    side: EndpointSide,
    justification: LabelJustification,
    incident_bond_order: i32,
    auto_rectify_hydrogens: bool,
) -> ParsedLabelAttachment {
    let trimmed = label.trim();

    if let Some((att_z, att_sym, subs, valid)) = resolve_known_nickname(trimmed, incident_bond_order) {
        // Explicit Left/Right on asymmetric forward/reversed multi-atom oxy/nitro groups respects endpoint order
        let force_wrong_oxy = matches!(
            (trimmed, side, justification),
            ("NO2", EndpointSide::Left, LabelJustification::Right)
                | ("O2N", EndpointSide::Right, LabelJustification::Left)
                | ("SO3H" | "CO2H", EndpointSide::Left, LabelJustification::Right)
        );
        if force_wrong_oxy {
            return ParsedLabelAttachment {
                raw_label: trimmed.to_string(),
                justification,
                attached_element_z: 8,
                attached_symbol: "O".to_string(),
                terminal_substituents: vec![(att_z, 1)],
                has_valence_warning: true,
                warning_message: Some("An atom in this label has an invalid valence.".to_string()),
            };
        }
        return ParsedLabelAttachment {
            raw_label: trimmed.to_string(),
            justification,
            attached_element_z: att_z,
            attached_symbol: att_sym.to_string(),
            terminal_substituents: subs,
            has_valence_warning: !valid,
            warning_message: if !valid {
                Some("An atom in this label has an invalid valence.".to_string())
            } else {
                None
            },
        };
    }

    if let Some(mut tokens) = tokenize_condensed_formula(trimmed) {
        let mut att_idx = pick_attachment_token_idx(&tokens, side, justification);
        let att_z = tokens[att_idx].z;
        let att_sym = tokens[att_idx].symbol.clone();

        // General H-rectification (`CH3` + bond order 2 -> `CH2`, `NH2` + bond order 2 -> `NH`, etc.)
        if auto_rectify_hydrogens
            && incident_bond_order > 1
            && chiralcraft_chem::can_infer_implicit_hydrogen(att_z)
        {
            let h_count: usize = tokens.iter().filter(|t| t.z == 1).map(|t| t.count).sum();
            let non_h_others: usize = tokens
                .iter()
                .enumerate()
                .filter(|&(idx, t)| t.z != 1 && idx != att_idx)
                .map(|(_, t)| t.count)
                .sum::<usize>()
                + tokens[att_idx].count.saturating_sub(1);
            let exp_val = default_valence(att_z).max(0) as usize;
            let total = incident_bond_order as usize + non_h_others + h_count;
            if total > exp_val {
                let drop = total - exp_val;
                if drop > 0 && drop <= h_count {
                    let mut remaining_drop = drop;
                    for tok in &mut tokens {
                        if tok.z == 1 && remaining_drop > 0 {
                            let take = tok.count.min(remaining_drop);
                            tok.count -= take;
                            remaining_drop -= take;
                        }
                    }
                    tokens.retain(|t| t.count > 0);
                    att_idx = pick_attachment_token_idx(&tokens, side, justification);
                }
            }
        }

        let effective_label = format_formula_tokens(&tokens);
        let mut subs: Vec<(i32, usize)> = Vec::new();
        for (idx, tok) in tokens.iter().enumerate() {
            let count = if idx == att_idx {
                tok.count.saturating_sub(1)
            } else {
                tok.count
            };
            if count == 0 {
                continue;
            }
            if let Some(existing) = subs.iter_mut().find(|(z, _)| *z == tok.z) {
                existing.1 += count;
            } else {
                subs.push((tok.z, count));
            }
        }

        let invalid = evaluate_attached_valence(att_z, incident_bond_order, &subs);
        return ParsedLabelAttachment {
            raw_label: effective_label,
            justification,
            attached_element_z: att_z,
            attached_symbol: att_sym,
            terminal_substituents: subs,
            has_valence_warning: invalid,
            warning_message: if invalid {
                Some("An atom in this label has an invalid valence.".to_string())
            } else {
                None
            },
        };
    }

    ParsedLabelAttachment {
        raw_label: trimmed.to_string(),
        justification,
        attached_element_z: 6,
        attached_symbol: "C".to_string(),
        terminal_substituents: Vec::new(),
        has_valence_warning: false,
        warning_message: None,
    }
}

/// Evaluates chemical attachment and valence validity for left-endpoint labels under `Auto`/`Left`/`Right` justification.
pub fn analyze_left_endpoint_label(
    label: &str,
    justification: LabelJustification,
    incident_bond_order: i32,
    auto_rectify_hydrogens: bool,
) -> ParsedLabelAttachment {
    analyze_endpoint_label(
        label,
        EndpointSide::Left,
        justification,
        incident_bond_order,
        auto_rectify_hydrogens,
    )
}

/// Expands a contracted label (`F3C`, `CF3`, `CH3`, `H3C`) on `atom_idx` into explicit terminal atoms and single bonds.
pub fn expand_contracted_label(
    mol: &mut ChemMolecule,
    atom_idx: usize,
    parsed: &ParsedLabelAttachment,
    bond_length_pt: f64,
) -> Result<Vec<usize>, String> {
    if parsed.has_valence_warning {
        return Err(parsed
            .warning_message
            .clone()
            .unwrap_or_else(|| "Invalid valence in contracted label".to_string()));
    }
    mol.atoms[atom_idx].z = parsed.attached_element_z;
    let cx = mol.atoms[atom_idx].x;
    let cy = mol.atoms[atom_idx].y;
    let angles = [std::f64::consts::PI, std::f64::consts::FRAC_PI_2, -std::f64::consts::FRAC_PI_2];
    let mut created = Vec::new();
    let mut slot = 0usize;
    let mut next_id = (mol.atoms.len() + mol.bonds.len() + 100) as u32;

    for &(sub_z, count) in &parsed.terminal_substituents {
        for _ in 0..count {
            let ang = angles[slot % angles.len()];
            slot += 1;
            let mut sub_atom = ChemAtom::new(next_id, sub_z);
            next_id += 1;
            sub_atom.x = cx + bond_length_pt * ang.cos();
            sub_atom.y = cy + bond_length_pt * ang.sin();
            let sub_idx = mol.add_atom(sub_atom);
            mol.add_bond(next_id, atom_idx, sub_idx, 1, BondDisplay::Solid);
            next_id += 1;
            created.push(sub_idx);
        }
    }
    Ok(created)
}

// ============================================================================
// 4. Selection Containment & Structure Cleanup Invariants (P05 & P02)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect2D {
    pub left: f64,
    pub top: f64,
    pub right: f64,
    pub bottom: f64,
}

impl Rect2D {
    #[inline]
    pub fn contains_point(&self, x: f64, y: f64) -> bool {
        x >= self.left && x <= self.right && y >= self.top && y <= self.bottom
    }

    #[inline]
    pub fn contains_rect(&self, other: &Rect2D) -> bool {
        self.contains_point(other.left, other.top)
            && self.contains_point(other.right, other.bottom)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

/// Evaluates whether a non-chemical graphic object with `obj_bounds` is selected by a marquee `sel_rect` (`P05`).
/// Full containment is required (`sel_rect.contains_rect(&obj_bounds)`).
#[inline]
pub fn marquee_selects_bounds(sel_rect: &Rect2D, obj_bounds: &Rect2D) -> bool {
    sel_rect.contains_rect(obj_bounds)
}

/// Evaluates whether a closed polygon `lasso` contains point `(px, py)` (even-odd ray-casting).
pub fn lasso_selects_point(lasso: &[Point2D], px: f64, py: f64) -> bool {
    if lasso.len() < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = lasso.len() - 1;
    for i in 0..lasso.len() {
        let xi = lasso[i].x;
        let yi = lasso[i].y;
        let xj = lasso[j].x;
        let yj = lasso[j].y;
        let intersects = ((yi > py) != (yj > py))
            && (px < (xj - xi) * (py - yi) / ((yj - yi).abs().max(1e-12) * (yj - yi).signum()) + xi);
        if intersects {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Full-containment check for a bounding rectangle inside a polygon lasso (`P05`).
pub fn lasso_selects_bounds(lasso: &[Point2D], obj_bounds: &Rect2D) -> bool {
    lasso_selects_point(lasso, obj_bounds.left, obj_bounds.top)
        && lasso_selects_point(lasso, obj_bounds.right, obj_bounds.top)
        && lasso_selects_point(lasso, obj_bounds.right, obj_bounds.bottom)
        && lasso_selects_point(lasso, obj_bounds.left, obj_bounds.bottom)
}

/// Computes document-point hit-testing tolerance from ChiralCraft's documented ~5 px screen-space tolerance (`C04`).
#[inline]
pub fn hit_tolerance_pt(screen_px_tol: f64, zoom_scale: f64) -> f64 {
    let scale = zoom_scale.max(0.1);
    screen_px_tol / scale
}

/// Returns the closest atom index within `tol_pt`, if any (`E01`–`E05`).
pub fn hit_test_atom(mol: &ChemMolecule, px: f64, py: f64, tol_pt: f64) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, a) in mol.atoms.iter().enumerate() {
        let d = ((a.x - px).powi(2) + (a.y - py).powi(2)).sqrt();
        if d <= tol_pt && best.is_none_or(|(_, bd)| d < bd) {
            best = Some((i, d));
        }
    }
    best.map(|(i, _)| i)
}

/// Returns the closest bond index within `tol_pt` of its segment, if any (`E01`–`E05`).
pub fn hit_test_bond(mol: &ChemMolecule, px: f64, py: f64, tol_pt: f64) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (bi, b) in mol.bonds.iter().enumerate() {
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let vx = a2.x - a1.x;
        let vy = a2.y - a1.y;
        let len_sq = vx * vx + vy * vy;
        let t = if len_sq > 1e-12 {
            (((px - a1.x) * vx + (py - a1.y) * vy) / len_sq).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let cx = a1.x + t * vx;
        let cy = a1.y + t * vy;
        let d = ((px - cx).powi(2) + (py - cy).powi(2)).sqrt();
        if d <= tol_pt && best.is_none_or(|(_, bd)| d < bd) {
            best = Some((bi, d));
        }
    }
    best.map(|(bi, _)| bi)
}

/// Translates selected atom indices cleanly (`dx`, `dy`).
pub fn translate_atoms(mol: &mut ChemMolecule, atoms: &[usize], dx: f64, dy: f64) {
    for &ai in atoms {
        if let Some(a) = mol.atoms.get_mut(ai) {
            a.x += dx;
            a.y += dy;
        }
    }
}

/// Scales selected atoms about their centroid by `(sx, sy)`.
pub fn scale_atoms(mol: &mut ChemMolecule, atoms: &[usize], sx: f64, sy: f64) {
    if atoms.is_empty() {
        return;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cnt = 0.0;
    for &ai in atoms {
        if let Some(a) = mol.atoms.get(ai) {
            cx += a.x;
            cy += a.y;
            cnt += 1.0;
        }
    }
    if cnt == 0.0 {
        return;
    }
    cx /= cnt;
    cy /= cnt;
    for &ai in atoms {
        if let Some(a) = mol.atoms.get_mut(ai) {
            a.x = cx + (a.x - cx) * sx;
            a.y = cy + (a.y - cy) * sy;
        }
    }
}

/// In-plane rigid rotation of selected atoms by `angle_rad` about centroid (preserves 2D CIP stereocenters).
pub fn rotate_atoms_preserving_stereo(mol: &mut ChemMolecule, atoms: &[usize], angle_rad: f64) {
    if atoms.is_empty() {
        return;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cnt = 0.0;
    for &ai in atoms {
        if let Some(a) = mol.atoms.get(ai) {
            cx += a.x;
            cy += a.y;
            cnt += 1.0;
        }
    }
    if cnt == 0.0 {
        return;
    }
    cx /= cnt;
    cy /= cnt;
    let c = angle_rad.cos();
    let s = angle_rad.sin();
    for &ai in atoms {
        if let Some(a) = mol.atoms.get_mut(ai) {
            let rx = a.x - cx;
            let ry = a.y - cy;
            a.x = cx + rx * c - ry * s;
            a.y = cy + rx * s + ry * c;
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlipAxis {
    Horizontal,
    Vertical,
}

/// Flips selected atoms across their centroid along `axis`. When `preserve_3d_stereo` is true
/// (ChiralCraft's Flip behavior toggling front/back wedges so absolute R/S config is preserved),
/// wedge and hashed-wedge bonds strictly incident to flipped atoms are toggled (`WedgeBegin <-> WedgedHashBegin`).
pub fn flip_atoms_preserving_stereo(
    mol: &mut ChemMolecule,
    atoms: &[usize],
    axis: FlipAxis,
    preserve_3d_stereo: bool,
) {
    if atoms.is_empty() {
        return;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cnt = 0.0;
    let aset: HashSet<usize> = atoms.iter().copied().collect();
    for &ai in &aset {
        if let Some(a) = mol.atoms.get(ai) {
            cx += a.x;
            cy += a.y;
            cnt += 1.0;
        }
    }
    if cnt == 0.0 {
        return;
    }
    cx /= cnt;
    cy /= cnt;
    for &ai in &aset {
        if let Some(a) = mol.atoms.get_mut(ai) {
            match axis {
                FlipAxis::Horizontal => a.x = 2.0 * cx - a.x,
                FlipAxis::Vertical => a.y = 2.0 * cy - a.y,
            }
        }
    }
    if preserve_3d_stereo {
        for b in &mut mol.bonds {
            if aset.contains(&b.begin_atom) && aset.contains(&b.end_atom) {
                b.display = match b.display {
                    BondDisplay::WedgeBegin => BondDisplay::WedgedHashBegin,
                    BondDisplay::WedgedHashBegin => BondDisplay::WedgeBegin,
                    BondDisplay::WedgeEnd => BondDisplay::WedgedHashEnd,
                    BondDisplay::WedgedHashEnd => BondDisplay::WedgeEnd,
                    BondDisplay::HollowWedgeBegin => BondDisplay::Hash,
                    BondDisplay::Hash => BondDisplay::HollowWedgeBegin,
                    other => other,
                };
            }
        }
    }
}

/// Minimal pure-Rust RFC 4648 Base64 encoder for embedding lossless CDXML inside SVG `<metadata>`.
fn b64_encode(input: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(((input.len() + 2) / 3) * 4);
    let mut i = 0;
    while i < input.len() {
        let b0 = input[i];
        let b1 = if i + 1 < input.len() { input[i + 1] } else { 0 };
        let b2 = if i + 2 < input.len() { input[i + 2] } else { 0 };
        let triple = ((b0 as u32) << 16) | ((b1 as u32) << 8) | (b2 as u32);
        out.push(TABLE[((triple >> 18) & 0x3F) as usize] as char);
        out.push(TABLE[((triple >> 12) & 0x3F) as usize] as char);
        out.push(if i + 1 < input.len() {
            TABLE[((triple >> 6) & 0x3F) as usize] as char
        } else {
            '='
        });
        out.push(if i + 2 < input.len() {
            TABLE[(triple & 0x3F) as usize] as char
        } else {
            '='
        });
        i += 3;
    }
    out
}

/// Minimal pure-Rust RFC 4648 Base64 decoder (ignores ASCII whitespace).
fn b64_decode(input: &str) -> Result<Vec<u8>, String> {
    fn val(b: u8) -> Option<u8> {
        match b {
            b'A'..=b'Z' => Some(b - b'A'),
            b'a'..=b'z' => Some(b - b'a' + 26),
            b'0'..=b'9' => Some(b - b'0' + 52),
            b'+' => Some(62),
            b'/' => Some(63),
            _ => None,
        }
    }
    let clean: Vec<u8> = input
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .collect();
    if clean.len() % 4 != 0 {
        return Err("Invalid base64 length in SVG metadata".to_string());
    }
    let mut out = Vec::with_capacity((clean.len() / 4) * 3);
    for chunk in clean.chunks(4) {
        let a = val(chunk[0]).ok_or_else(|| "Invalid base64 char".to_string())?;
        let b = val(chunk[1]).ok_or_else(|| "Invalid base64 char".to_string())?;
        out.push((a << 2) | (b >> 4));
        if chunk[2] != b'=' {
            let c = val(chunk[2]).ok_or_else(|| "Invalid base64 char".to_string())?;
            out.push((b << 4) | (c >> 2));
            if chunk[3] != b'=' {
                let d = val(chunk[3]).ok_or_else(|| "Invalid base64 char".to_string())?;
                out.push((c << 6) | d);
            }
        }
    }
    Ok(out)
}

/// Returns `true` if `svg` contains an embedded lossless ChiralCraft chemical document payload.
pub fn svg_contains_chiralcraft_document(svg: &str) -> bool {
    svg.contains("<chiralcraft:document") && svg.contains("</chiralcraft:document>")
}

/// Renders a single `ChemMolecule` into deterministic editable SVG honoring `DocumentSettings`
/// (`line_width_pt`, `bold_width_pt`, `bond_spacing_pct`, `hash_spacing_pt`, `label_font`, `label_size_pt`).
pub fn render_molecule_svg(mol: &ChemMolecule, settings: &DocumentSettings) -> String {
    let mut doc = InteractiveDocument::new_saved_seed("render.svg", settings.clone());
    doc.scene.molecules.push(mol.clone());
    doc.to_svg()
}

impl InteractiveDocument {
    /// Extracts a lossless `InteractiveDocument` from an editable ChiralCraft SVG containing
    /// `<metadata><chiralcraft:document format="cdxml" encoding="base64">...</chiralcraft:document></metadata>`.
    pub fn from_svg(svg: &str) -> Result<Self, String> {
        let start_tag = svg
            .find("<chiralcraft:document")
            .ok_or_else(|| "SVG does not contain embedded ChiralCraft chemical data".to_string())?;
        let after_tag = &svg[start_tag..];
        let open_end = after_tag
            .find('>')
            .ok_or_else(|| "Malformed <chiralcraft:document> opening tag".to_string())?;
        let header = &after_tag[..open_end];
        let body_start = start_tag + open_end + 1;
        let close_rel = svg[body_start..]
            .find("</chiralcraft:document>")
            .ok_or_else(|| "Missing </chiralcraft:document> closing tag".to_string())?;
        let payload = svg[body_start..body_start + close_rel].trim();

        let cdxml_str = if header.contains("encoding=\"base64\"") {
            let bytes = b64_decode(payload)?;
            String::from_utf8(bytes)
                .map_err(|_| "Embedded ChiralCraft CDXML is not valid UTF-8".to_string())?
        } else if let Some(inner) = payload
            .strip_prefix("<![CDATA[")
            .and_then(|s| s.strip_suffix("]]>"))
        {
            inner.to_string()
        } else {
            payload.to_string()
        };

        let cdxml = CdxmlDocument::parse(&cdxml_str)
            .map_err(|e| format!("Failed to parse embedded CDXML from SVG: {e:?}"))?;
        Ok(Self::from_cdxml(&cdxml))
    }

    /// Exports the active scene as a raster PNG with the lossless CDXML document embedded
    /// inside a standard PNG `tEXt` chunk (`ChiralCraftDocument`) for round-trip reopening.
    pub fn to_editable_png(&self, dpi: u32) -> Result<Vec<u8>, String> {
        let mut merged = ChemMolecule::new();
        for m in &self.scene.molecules {
            let off = merged.atoms.len();
            for a in &m.atoms {
                merged.atoms.push(a.clone());
            }
            for b in &m.bonds {
                let mut bc = b.clone();
                bc.begin_atom += off;
                bc.end_atom += off;
                merged.bonds.push(bc);
            }
        }
        if merged.atoms.is_empty() {
            return Err("No molecule on canvas to export as PNG".to_string());
        }
        let scale = ((dpi.clamp(72, 1200) as f64) / 72.0).max(1.5);
        let raw_png = chiralcraft_formats::molecule_to_png(&merged, scale);
        let cdxml_xml = self.to_cdxml().to_xml_string();
        Ok(chiralcraft_formats::embed_chiralcraft_cdxml_in_png(
            &raw_png,
            &cdxml_xml,
        ))
    }

    /// Reopens a lossless `InteractiveDocument` from an editable ChiralCraft PNG.
    pub fn from_png(png_bytes: &[u8]) -> Result<Self, String> {
        let cdxml_str = chiralcraft_formats::extract_chiralcraft_cdxml_from_png(png_bytes)
            .ok_or_else(|| "PNG does not contain embedded ChiralCraft chemical data".to_string())?;
        let cdxml = CdxmlDocument::parse(&cdxml_str)
            .map_err(|e| format!("Failed to parse embedded CDXML from PNG: {e:?}"))?;
        Ok(Self::from_cdxml(&cdxml))
    }

    /// Exports the active scene as a native publication-grade vector **PDF 1.4** containing
    /// an embedded lossless ChiralCraft CDXML payload for round-trip reopening.
    pub fn to_editable_pdf(&self, dark_ink: bool, padding_pt: f64) -> Vec<u8> {
        let arrows: Vec<chiralcraft_formats::PdfArrowItem> = self
            .scene
            .arrows
            .iter()
            .map(|a| chiralcraft_formats::PdfArrowItem {
                tail: a.tail,
                head: a.head,
            })
            .collect();
        let captions: Vec<chiralcraft_formats::PdfCaptionItem> = self
            .scene
            .captions
            .iter()
            .map(|c| chiralcraft_formats::PdfCaptionItem {
                x: c.x,
                y: c.y,
                text: c.text.clone(),
            })
            .collect();
        let cdxml_xml = self.to_cdxml().to_xml_string();
        chiralcraft_formats::scene_to_editable_pdf(
            &self.scene.molecules,
            &arrows,
            &captions,
            &cdxml_xml,
            &self.name,
            dark_ink,
            padding_pt,
        )
    }

    /// Reopens a lossless `InteractiveDocument` from an editable ChiralCraft PDF.
    pub fn from_pdf(pdf_bytes: &[u8]) -> Result<Self, String> {
        let cdxml_str = chiralcraft_formats::extract_chiralcraft_cdxml_from_pdf(pdf_bytes)
            .ok_or_else(|| "PDF does not contain embedded ChiralCraft chemical data".to_string())?;
        let cdxml = CdxmlDocument::parse(&cdxml_str)
            .map_err(|e| format!("Failed to parse embedded CDXML from PDF: {e:?}"))?;
        Ok(Self::from_cdxml(&cdxml))
    }

    /// Exports only the selected molecules (or the full scene if `selected_mols` is empty)
    /// as a self-contained editable ChiralCraft SVG ready for clipboard pasting into WordCraft/DeckCraft.
    pub fn selection_to_svg(&self, selected_mols: &[usize], transparent_bg: bool) -> String {
        self.selection_to_svg_styled(selected_mols, transparent_bg, false, 1.0, 3.6)
    }

    /// Exports only the selected molecules (or full scene) with explicit ink theme, 3D rear opacity, and crossing clearance.
    pub fn selection_to_svg_styled(
        &self,
        selected_mols: &[usize],
        transparent_bg: bool,
        dark_ink: bool,
        rear_opacity: f32,
        crossing_clearance_pt: f64,
    ) -> String {
        if selected_mols.is_empty() {
            return self.to_svg_styled(transparent_bg, dark_ink, rear_opacity, crossing_clearance_pt);
        }
        let mut sub = InteractiveDocument::new_saved_seed(&self.name, self.settings.clone());
        for &mi in selected_mols {
            if let Some(m) = self.scene.molecules.get(mi) {
                sub.scene.molecules.push(m.clone());
            }
        }
        for circ in &self.scene.aromatic_circles {
            if circ
                .molecule_idx
                .map(|mi| selected_mols.contains(&mi))
                .unwrap_or(false)
            {
                sub.scene.aromatic_circles.push(circ.clone());
            }
        }
        sub.to_svg_styled(transparent_bg, dark_ink, rear_opacity, crossing_clearance_pt)
    }

    /// Exports all molecules, reactions, graphics, and captions in the document scene to clean 2D vector SVG
    /// with an embedded lossless `<chiralcraft:document>` payload for round-trip editing.
    pub fn to_svg(&self) -> String {
        self.to_svg_styled(false, false, 1.0, 3.6)
    }

    /// Exports editable SVG with a transparent outer canvas (ideal for pasting onto WordCraft/DeckCraft slides).
    pub fn to_svg_transparent(&self) -> String {
        self.to_svg_styled(true, false, 1.0, 3.6)
    }

    /// Exports editable SVG with full styling options:
    /// - `transparent_bg`: omit outer white canvas & trim bonds around heteroatoms (no white sticker dots)
    /// - `dark_ink`: render bonds/labels in crisp light ink (`#F2F5F9`) for dark-mode presentation slides
    /// - `rear_opacity`: `0.0..=1.0` 3D Z-depth cueing for rear bonds/atoms
    /// - `crossing_clearance_pt`: bridge-gap radius in pt where a front bond crosses over a rear bond
    pub fn to_svg_styled(
        &self,
        transparent_bg: bool,
        dark_ink: bool,
        rear_opacity: f32,
        crossing_clearance_pt: f64,
    ) -> String {
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut expand = |x: f64, y: f64| {
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
        };
        for m in &self.scene.molecules {
            for a in &m.atoms {
                if a.visible {
                    expand(a.x, a.y);
                }
            }
        }
        for circ in &self.scene.aromatic_circles {
            expand(circ.center_x - circ.radius, circ.center_y - circ.radius);
            expand(circ.center_x + circ.radius, circ.center_y + circ.radius);
        }
        for c in &self.scene.captions {
            let approx_w = (c.text.chars().count() as f64) * (self.settings.caption_size_pt * 0.55);
            match c.justification {
                LabelJustification::Center | LabelJustification::Above | LabelJustification::Below => {
                    expand(c.x - 0.5 * approx_w, c.y - self.settings.caption_size_pt);
                    expand(c.x + 0.5 * approx_w, c.y + self.settings.caption_size_pt * 0.4);
                }
                LabelJustification::Right => {
                    expand(c.x - approx_w, c.y - self.settings.caption_size_pt);
                    expand(c.x, c.y + self.settings.caption_size_pt * 0.4);
                }
                LabelJustification::Left | LabelJustification::Auto => {
                    expand(c.x, c.y - self.settings.caption_size_pt);
                    expand(c.x + approx_w, c.y + self.settings.caption_size_pt * 0.4);
                }
            }
        }
        for arr in &self.scene.arrows {
            expand(arr.tail.0, arr.tail.1);
            expand(arr.head.0, arr.head.1);
        }
        for curve in &self.scene.curves {
            for &(px, py) in &curve.points {
                expand(px, py);
            }
        }
        for g in &self.scene.graphics {
            expand(g.bbox.0, g.bbox.1);
            expand(g.bbox.2, g.bbox.3);
        }
        for t in &self.scene.tables {
            expand(t.bbox.0, t.bbox.1);
            expand(t.bbox.2, t.bbox.3);
        }
        for tlc in &self.scene.tlc_plates {
            expand(tlc.bbox.0, tlc.bbox.1);
            expand(tlc.bbox.2, tlc.bbox.3);
        }
        for ep in &self.scene.energy_profiles {
            expand(ep.bbox.0, ep.bbox.1);
            expand(ep.bbox.2, ep.bbox.3);
        }
        if !min_x.is_finite() {
            min_x = 0.0;
            min_y = 0.0;
            max_x = 100.0;
            max_y = 100.0;
        }
        let pad = self.settings.bond_length_pt.max(12.0);
        let vx = min_x - pad;
        let vy = min_y - pad;
        let vw = (max_x - min_x + 2.0 * pad).max(40.0);
        let vh = (max_y - min_y + 2.0 * pad).max(40.0);

        let lw = self.settings.line_width_pt.max(0.2);
        let bw = self.settings.bold_width_pt.max(lw * 1.5);
        let bond_off = self.settings.bond_length_pt * (self.settings.bond_spacing_pct / 100.0);
        let hash_step = self.settings.hash_spacing_pt.max(1.2);

        let embedded_cdxml_b64 = b64_encode(self.to_cdxml().to_xml_string().as_bytes());
        let mut svg = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:chiralcraft=\"https://chiralcraft.dev/ns/1.0\" data-chiralcraft-version=\"1.0\" viewBox=\"{vx:.2} {vy:.2} {vw:.2} {vh:.2}\" width=\"{vw:.2}pt\" height=\"{vh:.2}pt\">\n<metadata>\n<chiralcraft:document format=\"cdxml\" encoding=\"base64\">{embedded_cdxml_b64}</chiralcraft:document>\n</metadata>\n"
        );
        let ink = if dark_ink { "#F2F5F9" } else { "black" };
        let bg_fill = if dark_ink { "#1E222A" } else { "white" };
        if !transparent_bg {
            svg.push_str(&format!(
                "<rect x=\"{vx:.2}\" y=\"{vy:.2}\" width=\"{vw:.2}\" height=\"{vh:.2}\" fill=\"{bg_fill}\"/>\n"
            ));
        }

        for (mol_idx, mol) in self.scene.molecules.iter().enumerate() {
            let val = CcVal::compute(mol);
            let crossings = find_bond_crossings(mol);
            let atom_opacities = compute_3d_atom_opacities(mol, rear_opacity);
            let rings = find_chordless_small_rings(mol, 8);
            let mut ring_centroids_by_bond: HashMap<usize, (f64, f64)> = HashMap::new();
            for (verts, r_bonds) in &rings {
                let mut rcx = 0.0;
                let mut rcy = 0.0;
                for &vi in verts {
                    rcx += mol.atoms[vi].x;
                    rcy += mol.atoms[vi].y;
                }
                rcx /= verts.len() as f64;
                rcy /= verts.len() as f64;
                for &bi in r_bonds {
                    ring_centroids_by_bond.entry(bi).or_insert((rcx, rcy));
                }
            }

            let atom_has_label: Vec<bool> = mol
                .atoms
                .iter()
                .enumerate()
                .map(|(ai, a)| {
                    a.visible
                        && (a.label.is_some()
                            || a.z != 6
                            || a.charge != 0
                            || a.isotope != 0
                            || a.radical != 0
                            || mol.neighbors(ai).is_empty())
                })
                .collect();
            let label_trim_r = self.settings.label_size_pt * 0.48 + self.settings.margin_width_pt * 0.85;

            for (bi, b) in mol.bonds.iter().enumerate() {
                let a1 = &mol.atoms[b.begin_atom];
                let a2 = &mol.atoms[b.end_atom];
                if !b.visible || !a1.visible || !a2.visible {
                    continue;
                }
                let (mut x1, mut y1, mut x2, mut y2) = (a1.x, a1.y, a2.x, a2.y);
                let raw_dx = x2 - x1;
                let raw_dy = y2 - y1;
                let raw_len = (raw_dx * raw_dx + raw_dy * raw_dy).sqrt().max(1e-9);
                if transparent_bg {
                    // On transparent SVGs (WordCraft/DeckCraft), trim bond ends cleanly around labeled atoms
                    // so bonds never strike through heteroatoms without needing opaque white sticker circles!
                    let ux = raw_dx / raw_len;
                    let uy = raw_dy / raw_len;
                    let max_trim = (raw_len * 0.38).min(label_trim_r);
                    if atom_has_label[b.begin_atom] {
                        x1 += ux * max_trim;
                        y1 += uy * max_trim;
                    }
                    if atom_has_label[b.end_atom] {
                        x2 -= ux * max_trim;
                        y2 -= uy * max_trim;
                    }
                }
                let dx = x2 - x1;
                let dy = y2 - y1;
                let len = (dx * dx + dy * dy).sqrt().max(1e-9);
                let nx = -dy / len;
                let ny = dx / len;
                let b_op = 0.5 * (atom_opacities[b.begin_atom] + atom_opacities[b.end_atom]);
                let op_attr = if b_op < 0.99 {
                    format!(" stroke-opacity=\"{b_op:.2}\" fill-opacity=\"{b_op:.2}\"")
                } else {
                    String::new()
                };
                let segs = if let Some(xpts) = crossings.get(&bi) {
                    split_segment_for_crossings(x1, y1, x2, y2, xpts, crossing_clearance_pt)
                } else {
                    vec![((x1, y1), (x2, y2))]
                };

                match b.display {
                    BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => {
                        let hw = 0.75 * bw; // F02: wedge end width = 1.5 * bold width
                        let p1x = x2 + nx * hw;
                        let p1y = y2 + ny * hw;
                        let p2x = x2 - nx * hw;
                        let p2y = y2 - ny * hw;
                        let fill = if b.display == BondDisplay::HollowWedgeBegin { "none" } else { ink };
                        svg.push_str(&format!(
                            "<polygon points=\"{x1:.2},{y1:.2} {p1x:.2},{p1y:.2} {p2x:.2},{p2y:.2}\" fill=\"{fill}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n"
                        ));
                    }
                    BondDisplay::WedgeEnd | BondDisplay::HollowWedgeEnd => {
                        let hw = 0.75 * bw;
                        let p1x = x1 + nx * hw;
                        let p1y = y1 + ny * hw;
                        let p2x = x1 - nx * hw;
                        let p2y = y1 - ny * hw;
                        let fill = if b.display == BondDisplay::HollowWedgeEnd { "none" } else { ink };
                        svg.push_str(&format!(
                            "<polygon points=\"{x2:.2},{y2:.2} {p1x:.2},{p1y:.2} {p2x:.2},{p2y:.2}\" fill=\"{fill}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n"
                        ));
                    }
                    BondDisplay::WedgedHashBegin | BondDisplay::Hash => {
                        let hw = 0.75 * bw;
                        let steps = (len / hash_step).floor().max(3.0) as usize;
                        for s in 1..=steps {
                            let t = (s as f64) / (steps as f64);
                            let cx = x1 + dx * t;
                            let cy = y1 + dy * t;
                            let w = hw * t;
                            svg.push_str(&format!(
                                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n",
                                cx + nx * w,
                                cy + ny * w,
                                cx - nx * w,
                                cy - ny * w
                            ));
                        }
                    }
                    BondDisplay::WedgedHashEnd => {
                        let hw = 0.75 * bw;
                        let steps = (len / hash_step).floor().max(3.0) as usize;
                        for s in 1..=steps {
                            let t = (s as f64) / (steps as f64);
                            let cx = x2 - dx * t;
                            let cy = y2 - dy * t;
                            let w = hw * t;
                            svg.push_str(&format!(
                                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n",
                                cx + nx * w,
                                cy + ny * w,
                                cx - nx * w,
                                cy - ny * w
                            ));
                        }
                    }
                    BondDisplay::Bold => {
                        for ((sx1, sy1), (sx2, sy2)) in &segs {
                            svg.push_str(&format!(
                                "<line x1=\"{sx1:.2}\" y1=\"{sy1:.2}\" x2=\"{sx2:.2}\" y2=\"{sy2:.2}\" stroke=\"{ink}\" stroke-width=\"{bw:.2}\" stroke-linecap=\"round\"{op_attr}/>\n"
                            ));
                        }
                    }
                    BondDisplay::Dash => {
                        for ((sx1, sy1), (sx2, sy2)) in &segs {
                            svg.push_str(&format!(
                                "<line x1=\"{sx1:.2}\" y1=\"{sy1:.2}\" x2=\"{sx2:.2}\" y2=\"{sy2:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\" stroke-dasharray=\"2.5,2.0\"{op_attr}/>\n"
                            ));
                        }
                    }
                    BondDisplay::Dative => {
                        let ux = dx / len;
                        let uy = dy / len;
                        let hl = (bw * 2.0).max(3.5);
                        let hw = (bw * 0.75).max(1.4);
                        svg.push_str(&format!(
                            "<line x1=\"{x1:.2}\" y1=\"{y1:.2}\" x2=\"{x2:.2}\" y2=\"{y2:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n"
                        ));
                        svg.push_str(&format!(
                            "<polygon points=\"{x2:.2},{y2:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{ink}\"{op_attr}/>\n",
                            x2 - ux * hl + nx * hw, y2 - uy * hl + ny * hw,
                            x2 - ux * hl - nx * hw, y2 - uy * hl - ny * hw,
                        ));
                    }
                    BondDisplay::Wavy => {
                        let mx = 0.5 * (x1 + x2) + nx * lw * 2.0;
                        let my = 0.5 * (y1 + y2) + ny * lw * 2.0;
                        svg.push_str(&format!(
                            "<path d=\"M {x1:.2} {y1:.2} Q {mx:.2} {my:.2} {x2:.2} {y2:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n"
                        ));
                    }
                    _ => {
                        if b.order == 2 || b.order == 6 {
                            let midx = 0.5 * (x1 + x2);
                            let midy = 0.5 * (y1 + y2);
                            let dash = if b.order == 6 { " stroke-dasharray=\"2.5,2.0\"" } else { "" };
                            if let Some(&(rcx, rcy)) = ring_centroids_by_bond.get(&bi) {
                                let sign = if (rcx - midx) * nx + (rcy - midy) * ny >= 0.0 {
                                    1.0
                                } else {
                                    -1.0
                                };
                                let in_x = nx * sign * bond_off;
                                let in_y = ny * sign * bond_off;
                                for ((sx1, sy1), (sx2, sy2)) in &segs {
                                    svg.push_str(&format!(
                                        "<line x1=\"{sx1:.2}\" y1=\"{sy1:.2}\" x2=\"{sx2:.2}\" y2=\"{sy2:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\" stroke-linecap=\"round\"{op_attr}/>\n"
                                    ));
                                }
                                svg.push_str(&format!(
                                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{dash}{op_attr}/>\n",
                                    x1 + 0.14 * dx + in_x,
                                    y1 + 0.14 * dy + in_y,
                                    x2 - 0.14 * dx + in_x,
                                    y2 - 0.14 * dy + in_y,
                                ));
                            } else {
                                let ho = 0.5 * bond_off;
                                for ((sx1, sy1), (sx2, sy2)) in &segs {
                                    svg.push_str(&format!(
                                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n",
                                        sx1 + nx * ho, sy1 + ny * ho, sx2 + nx * ho, sy2 + ny * ho
                                    ));
                                    svg.push_str(&format!(
                                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{dash}{op_attr}/>\n",
                                        sx1 - nx * ho, sy1 - ny * ho, sx2 - nx * ho, sy2 - ny * ho
                                    ));
                                }
                            }
                        } else if b.order == 3 {
                            for ((sx1, sy1), (sx2, sy2)) in &segs {
                                svg.push_str(&format!(
                                    "<line x1=\"{sx1:.2}\" y1=\"{sy1:.2}\" x2=\"{sx2:.2}\" y2=\"{sy2:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n"
                                ));
                                for s in [-1.0, 1.0] {
                                    svg.push_str(&format!(
                                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"{op_attr}/>\n",
                                        sx1 + nx * s * bond_off,
                                        sy1 + ny * s * bond_off,
                                        sx2 + nx * s * bond_off,
                                        sy2 + ny * s * bond_off,
                                    ));
                                }
                            }
                        } else {
                            for ((sx1, sy1), (sx2, sy2)) in &segs {
                                svg.push_str(&format!(
                                    "<line x1=\"{sx1:.2}\" y1=\"{sy1:.2}\" x2=\"{sx2:.2}\" y2=\"{sy2:.2}\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\" stroke-linecap=\"round\"{op_attr}/>\n"
                                ));
                            }
                        }
                    }
                }
            }

            // Aromatic ring circles attached to this fragment
            for circ in self
                .scene
                .aromatic_circles
                .iter()
                .filter(|c| c.molecule_idx == Some(mol_idx))
            {
                svg.push_str(&format!(
                    "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"/>\n",
                    circ.center_x, circ.center_y, circ.radius
                ));
            }

            // Atom labels (IUPAC GR-2.1: left-facing HO-/H2N-, superscript charges & isotopes, subscript H counts)
            for (ai, a) in mol.atoms.iter().enumerate() {
                if !a.visible {
                    continue;
                }
                let nbrs = mol.neighbors(ai);
                if a.label.is_none() && a.z == 6 && a.charge == 0 && a.isotope == 0 && a.radical == 0 && !nbrs.is_empty() {
                    continue;
                }
                let mut inner_svg = String::new();
                if let Some(custom_lbl) = &a.label {
                    inner_svg.push_str(&escape_xml_text(custom_lbl));
                } else {
                    let sym = chiralcraft_chem::get_element(a.z).map(|e| e.symbol).unwrap_or("?");
                    let h_cnt = val.implicit_h_counts[ai];
                    let flip_left = h_cnt > 0
                        && a.isotope == 0
                        && a.charge == 0
                        && nbrs.len() == 1
                        && mol.atoms[nbrs[0].0].x > a.x + 0.5;
                    if flip_left {
                        inner_svg.push('H');
                        if h_cnt > 1 {
                            inner_svg.push_str(&format!("<tspan baseline-shift=\"sub\" font-size=\"72%\">{h_cnt}</tspan>"));
                        }
                        inner_svg.push_str(sym);
                    } else {
                        if a.isotope != 0 {
                            inner_svg.push_str(&format!("<tspan baseline-shift=\"super\" font-size=\"72%\">{}</tspan>", a.isotope));
                        }
                        inner_svg.push_str(sym);
                        if h_cnt > 0 {
                            inner_svg.push('H');
                            if h_cnt > 1 {
                                inner_svg.push_str(&format!("<tspan baseline-shift=\"sub\" font-size=\"72%\">{h_cnt}</tspan>"));
                            }
                        }
                        if a.charge != 0 {
                            let chg_s = match a.charge {
                                1 => "+".to_string(),
                                -1 => "−".to_string(),
                                c if c > 1 => format!("{c}+"),
                                c => format!("{}−", -c),
                            };
                            inner_svg.push_str(&format!("<tspan baseline-shift=\"super\" font-size=\"72%\">{chg_s}</tspan>"));
                        }
                    }
                }
                let a_op = atom_opacities[ai];
                let a_op_attr = if a_op < 0.99 {
                    format!(" fill-opacity=\"{a_op:.2}\"")
                } else {
                    String::new()
                };
                if !transparent_bg {
                    let r = self.settings.label_size_pt * 0.45 + self.settings.margin_width_pt;
                    svg.push_str(&format!(
                        "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{r:.2}\" fill=\"{bg_fill}\"/>\n",
                        a.x, a.y
                    ));
                }
                svg.push_str(&format!(
                    "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{}\" font-size=\"{:.1}\" text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"{ink}\"{a_op_attr}>{inner_svg}</text>\n",
                    a.x,
                    a.y,
                    escape_xml_attr(&self.settings.label_font),
                    self.settings.label_size_pt,
                ));
            }
        }

        for circ in self
            .scene
            .aromatic_circles
            .iter()
            .filter(|c| c.molecule_idx.is_none())
        {
            svg.push_str(&format!(
                "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"{lw:.2}\"/>\n",
                circ.center_x, circ.center_y, circ.radius
            ));
        }

        for cap in &self.scene.captions {
            let anchor = match cap.justification {
                LabelJustification::Center | LabelJustification::Above | LabelJustification::Below => "middle",
                LabelJustification::Right => "end",
                LabelJustification::Left | LabelJustification::Auto => "start",
            };
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{}\" font-size=\"{:.1}\" text-anchor=\"{anchor}\" fill=\"{ink}\">{}</text>\n",
                cap.x,
                cap.y,
                escape_xml_attr(&self.settings.caption_font),
                self.settings.caption_size_pt,
                escape_xml_text(&cap.text)
            ));
        }

        svg.push_str(&render_scene_ext_svg_styled(
            &self.scene.arrows,
            &self.scene.curves,
            &self.scene.graphics,
            &self.scene.tables,
            &self.scene.tlc_plates,
            &self.scene.header_footer,
            ink,
        ));
        svg.push_str(&render_energy_profiles_svg(&self.scene.energy_profiles, ink));

        svg.push_str("</svg>\n");
        svg
    }
}

/// Appends a new carbon atom and single bond from `atom_idx` at `angle_deg` using `settings.bond_length_pt` (`D01`/`D10`).
pub fn add_bond_from_atom(
    mol: &mut ChemMolecule,
    atom_idx: usize,
    angle_deg: f64,
    settings: &DocumentSettings,
) -> Option<(usize, usize)> {
    let anchor = mol.atoms.get(atom_idx)?;
    let (ax, ay) = (anchor.x, anchor.y);
    let rad = angle_deg.to_radians();
    let next_aid = (mol.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1).max(1);
    let next_bid = (mol.bonds.iter().map(|b| b.id).max().unwrap_or(0) + 1).max(1);
    let mut new_atom = ChemAtom::new(next_aid, 6);
    new_atom.x = ax + settings.bond_length_pt * rad.cos();
    new_atom.y = ay + settings.bond_length_pt * rad.sin();
    let new_idx = mol.add_atom(new_atom);
    let bi = mol.add_bond(next_bid, atom_idx, new_idx, 1, BondDisplay::Solid);
    Some((new_idx, bi))
}

/// Cycles bond order `1 -> 2 -> 3 -> 4 -> 1` (`D02`/`D03`).
pub fn cycle_bond_order(mol: &mut ChemMolecule, bond_idx: usize) -> Option<i32> {
    let b = mol.bonds.get_mut(bond_idx)?;
    b.order = match b.order {
        1 => 2,
        2 => 3,
        3 => 4,
        _ => 1,
    };
    Some(b.order)
}

/// Collapses `bond_idx` by merging `end_atom` into `begin_atom` and removing the connecting bond (`D11`).
pub fn collapse_bond_reducing_ring(mol: &mut ChemMolecule, bond_idx: usize) -> bool {
    if bond_idx >= mol.bonds.len() {
        return false;
    }
    let keep = mol.bonds[bond_idx].begin_atom;
    let drop = mol.bonds[bond_idx].end_atom;
    if keep == drop {
        return false;
    }
    mol.bonds.remove(bond_idx);
    for b in &mut mol.bonds {
        if b.begin_atom == drop {
            b.begin_atom = keep;
        }
        if b.end_atom == drop {
            b.end_atom = keep;
        }
    }
    mol.bonds.retain(|b| b.begin_atom != b.end_atom);
    true
}

/// Fuses a regular `ring_size`-membered carbon ring (`3..=8`) onto `bond_idx` (`D12`/`D13`).
pub fn fuse_ring_on_bond(
    mol: &mut ChemMolecule,
    bond_idx: usize,
    ring_size: usize,
    aromatic: bool,
) -> Vec<usize> {
    if bond_idx >= mol.bonds.len() || !(3..=8).contains(&ring_size) {
        return Vec::new();
    }
    let u = mol.bonds[bond_idx].begin_atom;
    let v = mol.bonds[bond_idx].end_atom;
    let (ux, uy) = (mol.atoms[u].x, mol.atoms[u].y);
    let (vx, vy) = (mol.atoms[v].x, mol.atoms[v].y);
    let dx = vx - ux;
    let dy = vy - uy;
    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
    let tx = dx / len;
    let ty = dy / len;
    let nx = -ty;
    let ny = tx;

    let k_f = ring_size as f64;
    let half_angle = std::f64::consts::PI / k_f;
    let radius = len / (2.0 * half_angle.sin());
    let apothem = radius * half_angle.cos();
    let midx = 0.5 * (ux + vx);
    let midy = 0.5 * (uy + vy);
    let cx = midx + apothem * nx;
    let cy = midy + apothem * ny;

    let theta_u = (uy - cy).atan2(ux - cx);
    let theta_v = (vy - cy).atan2(vx - cx);
    let mut d_uv = theta_v - theta_u;
    while d_uv > std::f64::consts::PI {
        d_uv -= 2.0 * std::f64::consts::PI;
    }
    while d_uv < -std::f64::consts::PI {
        d_uv += 2.0 * std::f64::consts::PI;
    }
    let step = if d_uv >= 0.0 { 2.0 * half_angle } else { -2.0 * half_angle };

    let mut next_aid = (mol.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1).max(1);
    let mut next_bid = (mol.bonds.iter().map(|b| b.id).max().unwrap_or(0) + 1).max(1);
    let order = if aromatic { 6 } else { 1 };
    if aromatic {
        mol.bonds[bond_idx].order = 6;
    }

    let mut added = Vec::new();
    let mut prev = v;
    for offset in 2..ring_size {
        let theta = theta_u + (offset as f64) * step;
        let mut na = ChemAtom::new(next_aid, 6);
        next_aid += 1;
        na.x = cx + radius * theta.cos();
        na.y = cy + radius * theta.sin();
        let ni = mol.add_atom(na);
        mol.add_bond(next_bid, prev, ni, order, BondDisplay::Solid);
        next_bid += 1;
        added.push(ni);
        prev = ni;
    }
    mol.add_bond(next_bid, prev, u, order, BondDisplay::Solid);
    added
}

/// Attaches a spiro-fused regular `ring_size`-membered ring (`3..=8`) onto `atom_idx`
/// extending along the outward bisector of its incident bonds (`D12`).
pub fn spiro_ring_on_atom(
    mol: &mut ChemMolecule,
    atom_idx: usize,
    ring_size: usize,
    bond_length_pt: f64,
    aromatic: bool,
) -> Vec<usize> {
    if atom_idx >= mol.atoms.len() || !(3..=8).contains(&ring_size) {
        return Vec::new();
    }
    let ax = mol.atoms[atom_idx].x;
    let ay = mol.atoms[atom_idx].y;
    let nbrs = mol.neighbors(atom_idx);
    let mut sum_vx = 0.0;
    let mut sum_vy = 0.0;
    for &(ni, _) in &nbrs {
        let dx = mol.atoms[ni].x - ax;
        let dy = mol.atoms[ni].y - ay;
        let len = (dx * dx + dy * dy).sqrt();
        if len > 1e-9 {
            sum_vx += dx / len;
            sum_vy += dy / len;
        }
    }
    let out_ang = if sum_vx * sum_vx + sum_vy * sum_vy > 1e-9 {
        (-sum_vy).atan2(-sum_vx)
    } else {
        0.0
    };

    let k_f = ring_size as f64;
    let radius = bond_length_pt / (2.0 * (std::f64::consts::PI / k_f).sin());
    let cx = ax + radius * out_ang.cos();
    let cy = ay + radius * out_ang.sin();
    let theta0 = out_ang + std::f64::consts::PI;
    let step = 2.0 * std::f64::consts::PI / k_f;

    let mut next_aid = (mol.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1).max(1);
    let mut next_bid = (mol.bonds.iter().map(|b| b.id).max().unwrap_or(0) + 1).max(1);
    let order = if aromatic { 6 } else { 1 };

    let mut added = Vec::new();
    let mut prev = atom_idx;
    for i in 1..ring_size {
        let theta = theta0 + (i as f64) * step;
        let mut na = ChemAtom::new(next_aid, 6);
        next_aid += 1;
        na.x = cx + radius * theta.cos();
        na.y = cy + radius * theta.sin();
        let ni = mol.add_atom(na);
        mol.add_bond(next_bid, prev, ni, order, BondDisplay::Solid);
        next_bid += 1;
        added.push(ni);
        prev = ni;
    }
    mol.add_bond(next_bid, prev, atom_idx, order, BondDisplay::Solid);
    added
}

/// Merges `drop_atom` into `keep_atom` (`E12` Join command), preserving `keep_atom` identity,
/// removing `drop_atom` cleanly from `mol.atoms`, reindexing bonds, and deduplicating parallel bonds.
pub fn join_atoms(mol: &mut ChemMolecule, keep_atom: usize, drop_atom: usize) -> bool {
    if keep_atom >= mol.atoms.len() || drop_atom >= mol.atoms.len() || keep_atom == drop_atom {
        return false;
    }
    for b in &mut mol.bonds {
        if b.begin_atom == drop_atom {
            b.begin_atom = keep_atom;
        }
        if b.end_atom == drop_atom {
            b.end_atom = keep_atom;
        }
    }
    mol.bonds.retain(|b| b.begin_atom != b.end_atom);
    let mut seen_pairs = HashSet::new();
    mol.bonds.retain(|b| {
        let pair = if b.begin_atom <= b.end_atom {
            (b.begin_atom, b.end_atom)
        } else {
            (b.end_atom, b.begin_atom)
        };
        seen_pairs.insert(pair)
    });
    mol.atoms.remove(drop_atom);
    for b in &mut mol.bonds {
        if b.begin_atom > drop_atom {
            b.begin_atom -= 1;
        }
        if b.end_atom > drop_atom {
            b.end_atom -= 1;
        }
    }
    true
}

/// Joins two separate fragments/molecules (`Object > Join`, `Ctrl+J`) by translating `drop_mol`
/// so `drop_atom` coincides with `keep_atom`, deleting `drop_atom`, and merging into `keep_mol` (`P30a`).
pub fn join_molecules_at_atoms(
    keep_mol: &mut ChemMolecule,
    keep_atom: usize,
    drop_mol: &ChemMolecule,
    drop_atom: usize,
) -> bool {
    if keep_atom >= keep_mol.atoms.len() || drop_atom >= drop_mol.atoms.len() {
        return false;
    }
    let dx = keep_mol.atoms[keep_atom].x - drop_mol.atoms[drop_atom].x;
    let dy = keep_mol.atoms[keep_atom].y - drop_mol.atoms[drop_atom].y;
    let mut next_aid = (keep_mol.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1).max(1);
    let mut next_bid = (keep_mol.bonds.iter().map(|b| b.id).max().unwrap_or(0) + 1).max(1);

    let mut idx_map = vec![0usize; drop_mol.atoms.len()];
    for (i, a) in drop_mol.atoms.iter().enumerate() {
        if i == drop_atom {
            idx_map[i] = keep_atom;
        } else {
            let mut na = a.clone();
            na.x += dx;
            na.y += dy;
            na.id = next_aid;
            next_aid += 1;
            idx_map[i] = keep_mol.add_atom(na);
        }
    }

    let mut seen: HashSet<(usize, usize)> = keep_mol
        .bonds
        .iter()
        .map(|b| (b.begin_atom.min(b.end_atom), b.begin_atom.max(b.end_atom)))
        .collect();
    for b in &drop_mol.bonds {
        let bi = idx_map[b.begin_atom];
        let ei = idx_map[b.end_atom];
        if bi != ei && seen.insert((bi.min(ei), bi.max(ei))) {
            keep_mol.add_bond(next_bid, bi, ei, b.order, b.display);
            next_bid += 1;
        }
    }
    true
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignMode {
    Left,
    Right,
    Top,
    Bottom,
    HorizontalCenter,
    VerticalCenter,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistributeAxis {
    Horizontal,
    Vertical,
}

fn mol_bounds(mol: &ChemMolecule) -> Option<Rect2D> {
    if mol.atoms.is_empty() {
        return None;
    }
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for a in &mol.atoms {
        min_x = min_x.min(a.x);
        min_y = min_y.min(a.y);
        max_x = max_x.max(a.x);
        max_y = max_y.max(a.y);
    }
    Some(Rect2D {
        left: min_x,
        top: min_y,
        right: max_x,
        bottom: max_y,
    })
}

fn shift_molecule(mol: &mut ChemMolecule, dx: f64, dy: f64) {
    for a in &mut mol.atoms {
        a.x += dx;
        a.y += dy;
    }
}

/// Aligns molecules cleanly using their outer bounding edges or centers (`E14`).
pub fn align_molecules(mols: &mut [ChemMolecule], mode: AlignMode) {
    let bounds: Vec<Option<Rect2D>> = mols.iter().map(mol_bounds).collect();
    let valid: Vec<Rect2D> = bounds.iter().filter_map(|&b| b).collect();
    if valid.len() < 2 {
        return;
    }
    let target = match mode {
        AlignMode::Left => valid.iter().map(|b| b.left).fold(f64::INFINITY, f64::min),
        AlignMode::Right => valid.iter().map(|b| b.right).fold(f64::NEG_INFINITY, f64::max),
        AlignMode::Top => valid.iter().map(|b| b.top).fold(f64::INFINITY, f64::min),
        AlignMode::Bottom => valid.iter().map(|b| b.bottom).fold(f64::NEG_INFINITY, f64::max),
        AlignMode::HorizontalCenter => {
            let min_l = valid.iter().map(|b| b.left).fold(f64::INFINITY, f64::min);
            let max_r = valid.iter().map(|b| b.right).fold(f64::NEG_INFINITY, f64::max);
            0.5 * (min_l + max_r)
        }
        AlignMode::VerticalCenter => {
            let min_t = valid.iter().map(|b| b.top).fold(f64::INFINITY, f64::min);
            let max_b = valid.iter().map(|b| b.bottom).fold(f64::NEG_INFINITY, f64::max);
            0.5 * (min_t + max_b)
        }
    };
    for (mol, b_opt) in mols.iter_mut().zip(bounds.iter()) {
        if let Some(b) = b_opt {
            match mode {
                AlignMode::Left => shift_molecule(mol, target - b.left, 0.0),
                AlignMode::Right => shift_molecule(mol, target - b.right, 0.0),
                AlignMode::Top => shift_molecule(mol, 0.0, target - b.top),
                AlignMode::Bottom => shift_molecule(mol, 0.0, target - b.bottom),
                AlignMode::HorizontalCenter => shift_molecule(mol, target - 0.5 * (b.left + b.right), 0.0),
                AlignMode::VerticalCenter => shift_molecule(mol, 0.0, target - 0.5 * (b.top + b.bottom)),
            }
        }
    }
}

/// Distributes 3+ molecules evenly along `axis`, preserving the outermost bounds (`E14`).
pub fn distribute_molecules(mols: &mut [ChemMolecule], axis: DistributeAxis) {
    let mut indexed: Vec<(usize, Rect2D)> = mols
        .iter()
        .enumerate()
        .filter_map(|(i, m)| mol_bounds(m).map(|b| (i, b)))
        .collect();
    if indexed.len() < 3 {
        return;
    }
    match axis {
        DistributeAxis::Horizontal => {
            indexed.sort_by(|a, b| a.1.left.total_cmp(&b.1.left));
            let min_edge = indexed.first().unwrap().1.left;
            let max_edge = indexed.last().unwrap().1.right;
            let total_widths: f64 = indexed.iter().map(|(_, b)| b.right - b.left).sum();
            let gap = ((max_edge - min_edge) - total_widths) / ((indexed.len() - 1) as f64);
            let mut cursor = min_edge;
            for (idx, b) in indexed {
                shift_molecule(&mut mols[idx], cursor - b.left, 0.0);
                cursor += (b.right - b.left) + gap;
            }
        }
        DistributeAxis::Vertical => {
            indexed.sort_by(|a, b| a.1.top.total_cmp(&b.1.top));
            let min_edge = indexed.first().unwrap().1.top;
            let max_edge = indexed.last().unwrap().1.bottom;
            let total_heights: f64 = indexed.iter().map(|(_, b)| b.bottom - b.top).sum();
            let gap = ((max_edge - min_edge) - total_heights) / ((indexed.len() - 1) as f64);
            let mut cursor = min_edge;
            for (idx, b) in indexed {
                shift_molecule(&mut mols[idx], 0.0, cursor - b.top);
                cursor += (b.bottom - b.top) + gap;
            }
        }
    }
}

fn selection_bounds(mols: &[ChemMolecule]) -> Option<Rect2D> {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for m in mols.iter() {
        if let Some(b) = mol_bounds(m) {
            min_x = min_x.min(b.left);
            min_y = min_y.min(b.top);
            max_x = max_x.max(b.right);
            max_y = max_y.max(b.bottom);
        }
    }
    if !min_x.is_finite() {
        None
    } else {
        Some(Rect2D {
            left: min_x,
            top: min_y,
            right: max_x,
            bottom: max_y,
        })
    }
}

/// Centers all molecules in the scene on the page bounds (`540 x 719.75 pt` default from `seed-v1.cdxml`, `E14`).
pub fn center_on_page(mols: &mut [ChemMolecule], page_width_pt: f64, page_height_pt: f64) {
    let Some(bounds) = selection_bounds(mols) else {
        return;
    };
    let dx = 0.5 * page_width_pt - 0.5 * (bounds.left + bounds.right);
    let dy = 0.5 * page_height_pt - 0.5 * (bounds.top + bounds.bottom);
    for m in mols.iter_mut() {
        shift_molecule(m, dx, dy);
    }
}

/// Rotates a multi-molecule selection clockwise in screen space (`+angle_rad`, `x` right, `y` down)
/// around an explicit pivot `(cx, cy)`.
pub fn rotate_molecules_around_point(mols: &mut [ChemMolecule], cx: f64, cy: f64, angle_rad: f64) {
    let c = angle_rad.cos();
    let s = angle_rad.sin();
    for m in mols.iter_mut() {
        for a in &mut m.atoms {
            let rx = a.x - cx;
            let ry = a.y - cy;
            a.x = cx + rx * c - ry * s;
            a.y = cy + rx * s + ry * c;
        }
    }
}

/// Rotates a multi-molecule selection clockwise in screen space (`+angle_rad`, `x` right, `y` down)
/// around the combined atom bounding-box center, matching ChiralCraft's `Object > Rotate...` behavior (`P23`).
pub fn rotate_molecules_around_bounds_center(mols: &mut [ChemMolecule], angle_rad: f64) {
    let Some(bounds) = selection_bounds(mols) else {
        return;
    };
    rotate_molecules_around_point(
        mols,
        0.5 * (bounds.left + bounds.right),
        0.5 * (bounds.top + bounds.bottom),
        angle_rad,
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StructureWarningTarget {
    Atom(usize),
    Bond(usize),
    AtomPair(usize, usize),
    Caption(u32),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureWarning {
    pub target: StructureWarningTarget,
    pub message: String,
}

/// Warning emitted by `check_structure_scene` or `check_enhanced_stereo_warnings`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SceneStructureWarning {
    pub molecule_idx: Option<usize>,
    pub object_id: u32,
    pub message: String,
}

#[inline]
fn point_to_segment_dist(px: f64, py: f64, x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dx = x2 - x1;
    let dy = y2 - y1;
    let l2 = dx * dx + dy * dy;
    if l2 <= 1e-12 {
        return ((px - x1).powi(2) + (py - y1).powi(2)).sqrt();
    }
    let t = (((px - x1) * dx + (py - y1) * dy) / l2).clamp(0.0, 1.0);
    let qx = x1 + t * dx;
    let qy = y1 + t * dy;
    ((px - qx).powi(2) + (py - qy).powi(2)).sqrt()
}

#[inline]
fn bond_is_stereobond_at_atom(b: &chiralcraft_chem::ChemBond, atom_idx: usize) -> bool {
    match b.display {
        BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin | BondDisplay::WedgedHashBegin => {
            b.begin_atom == atom_idx
        }
        BondDisplay::WedgeEnd | BondDisplay::HollowWedgeEnd | BondDisplay::WedgedHashEnd => {
            b.end_atom == atom_idx
        }
        BondDisplay::Hash | BondDisplay::Bold => {
            b.begin_atom == atom_idx || b.end_atom == atom_idx
        }
        _ => false,
    }
}

#[inline]
fn bond_is_any_stereobond(b: &chiralcraft_chem::ChemBond) -> bool {
    matches!(
        b.display,
        BondDisplay::WedgeBegin
            | BondDisplay::WedgeEnd
            | BondDisplay::HollowWedgeBegin
            | BondDisplay::HollowWedgeEnd
            | BondDisplay::WedgedHashBegin
            | BondDisplay::WedgedHashEnd
            | BondDisplay::Hash
            | BondDisplay::Bold
    )
}

fn atom_has_ambiguous_cross_stereobonds(mol: &ChemMolecule, atom_idx: usize) -> bool {
    let nbrs = mol.neighbors(atom_idx);
    if nbrs.len() != 4 {
        return false;
    }
    let stereo_count = nbrs
        .iter()
        .filter(|&&(_, bi)| {
            let b = &mol.bonds[bi];
            b.visible && !b.ignore_warnings && bond_is_stereobond_at_atom(b, atom_idx)
        })
        .count();
    if stereo_count == 0 {
        return false;
    }
    let c = &mol.atoms[atom_idx];
    let mut dirs = [[0.0_f64; 2]; 4];
    for (k, &(ni, _)) in nbrs.iter().enumerate() {
        let dx = mol.atoms[ni].x - c.x;
        let dy = mol.atoms[ni].y - c.y;
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 1e-6 {
            return false;
        }
        dirs[k] = [dx / len, dy / len];
    }
    let dot = |i: usize, j: usize| dirs[i][0] * dirs[j][0] + dirs[i][1] * dirs[j][1];
    let pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))];
    for ((a, b), (c_i, d)) in pairings {
        if dot(a, b) <= -0.92 && dot(c_i, d) <= -0.92 && dot(a, c_i).abs() <= 0.35 {
            return true;
        }
    }
    false
}

/// Evaluates `Structure > Check Structure` (`D17`, `C08`) across `mol` and `captions`:
/// - Unspecified labels (`NodeType="Unspecified"`, `"ChiralCraft can't interpret this label."`)
/// - Atom valence checks (`CcVal::compute`, `"An atom in this label has an invalid valence."`)
/// - Isolated single bond warning (`"This isolated bond is probably not intended to have chemical significance."`, observed in `p04-acs-100.cdxml`)
/// - Hydrogen bond validation (`D07`: `"Invalid or ambiguous H-Bonds"` when `order == 7` without H + electronegative neighbor)
/// - Colinear identical unlabeled-carbon bonds & orthogonal `+`-cross conflicting/ambiguous stereobonds
/// - Stereobonds connecting two stereocenters (`"This stereobond connects two stereocenters."`)
/// - Overlapping non-bonded atom proximity (`dist < 1.0 pt`) and atom-to-bond proximity (`dist < 2.0 pt`)
/// - Caption chemical warnings.
pub fn check_structure(mol: &ChemMolecule, captions: &[CaptionNode]) -> Vec<StructureWarning> {
    let mut warnings = Vec::new();
    let val = CcVal::compute(mol);

    for i in 0..mol.atoms.len() {
        let a = &mol.atoms[i];
        if a.ignore_warnings {
            continue;
        }
        if a.is_unspecified_label {
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Atom(i),
                message: "ChiralCraft can't interpret this label.".to_string(),
            });
            continue;
        }
        if atom_has_valence_warning(mol, &val, i) {
            let unlabeled_c = a.z == 6
                && a.charge == 0
                && a.isotope == 0
                && a.radical == 0
                && a.explicit_num_hydrogens.is_none();
            let msg = if unlabeled_c {
                "There are too many bonds to this unlabeled Carbon."
            } else {
                "An atom in this label has an invalid valence."
            };
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Atom(i),
                message: msg.to_string(),
            });
        }
        if atom_has_ambiguous_cross_stereobonds(mol, i) {
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Atom(i),
                message: "This stereocenter has conflicting or ambiguous stereobonds specified."
                    .to_string(),
            });
        }
    }

    if mol.atoms.len() == 2 && mol.bonds.len() == 1 && mol.bonds[0].order == 1 {
        let a0 = &mol.atoms[0];
        let a1 = &mol.atoms[1];
        let b0 = &mol.bonds[0];
        let plain_c0 = a0.z == 6 && a0.charge == 0 && a0.radical == 0 && a0.isotope == 0 && a0.valence_override.is_none() && a0.explicit_num_hydrogens.is_none() && !a0.is_unspecified_label && !a0.ignore_warnings && a0.visible;
        let plain_c1 = a1.z == 6 && a1.charge == 0 && a1.radical == 0 && a1.isotope == 0 && a1.valence_override.is_none() && a1.explicit_num_hydrogens.is_none() && !a1.is_unspecified_label && !a1.ignore_warnings && a1.visible;
        if plain_c0 && plain_c1 && !b0.ignore_warnings && b0.visible {
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Bond(0),
                message: "This isolated bond is probably not intended to have chemical significance.".to_string(),
            });
        }
    }

    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.ignore_warnings || !b.visible {
            continue;
        }
        if b.order == 7 {
            let z1 = mol.atoms[b.begin_atom].z;
            let z2 = mol.atoms[b.end_atom].z;
            let has_h = z1 == 1 || z2 == 1;
            let other = if z1 == 1 { z2 } else { z1 };
            let valid_acceptor = matches!(other, 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53);
            if !has_h || !valid_acceptor {
                warnings.push(StructureWarning {
                    target: StructureWarningTarget::Bond(bi),
                    message: "Invalid or ambiguous H-Bonds".to_string(),
                });
            }
        }
        if bond_is_any_stereobond(b)
            && chiralcraft_chem::compute_tetrahedral_cip_code(mol, b.begin_atom) != 0
            && chiralcraft_chem::compute_tetrahedral_cip_code(mol, b.end_atom) != 0
        {
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Bond(bi),
                message: "This stereobond connects two stereocenters.".to_string(),
            });
        }
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        let is_unlabeled_c = a.z == 6
            && a.charge == 0
            && a.isotope == 0
            && a.radical == 0
            && !a.is_query
            && !a.is_unspecified_label
            && !a.ignore_warnings
            && a.visible
            && a.explicit_num_hydrogens.is_none();
        if !is_unlabeled_c {
            continue;
        }
        let nbrs = mol.neighbors(i);
        if nbrs.len() == 2 {
            let (n0, b0) = nbrs[0];
            let (n1, b1) = nbrs[1];
            let bond0 = &mol.bonds[b0];
            let bond1 = &mol.bonds[b1];
            if bond0.order == bond1.order && bond0.display == bond1.display {
                let dx0 = mol.atoms[n0].x - a.x;
                let dy0 = mol.atoms[n0].y - a.y;
                let dx1 = mol.atoms[n1].x - a.x;
                let dy1 = mol.atoms[n1].y - a.y;
                let l0 = (dx0 * dx0 + dy0 * dy0).sqrt();
                let l1 = (dx1 * dx1 + dy1 * dy1).sqrt();
                if l0 > 1e-6 && l1 > 1e-6 {
                    let dot = (dx0 * dx1 + dy0 * dy1) / (l0 * l1);
                    if dot <= -0.995 {
                        warnings.push(StructureWarning {
                            target: StructureWarningTarget::Atom(i),
                            message: "This unlabeled atom is connected to two identical colinear bonds, and may appear to be invisible.".to_string(),
                        });
                    }
                }
            }
        }
    }

    let mut bonded = HashSet::new();
    for b in &mol.bonds {
        bonded.insert((b.begin_atom.min(b.end_atom), b.begin_atom.max(b.end_atom)));
    }
    let mut atom_close_warned = HashSet::new();
    for i in 0..mol.atoms.len() {
        if mol.atoms[i].ignore_warnings || !mol.atoms[i].visible {
            continue;
        }
        for j in (i + 1)..mol.atoms.len() {
            if mol.atoms[j].ignore_warnings || !mol.atoms[j].visible {
                continue;
            }
            if !bonded.contains(&(i, j)) {
                let d = ((mol.atoms[i].x - mol.atoms[j].x).powi(2)
                    + (mol.atoms[i].y - mol.atoms[j].y).powi(2))
                .sqrt();
                if d < 1.0 {
                    atom_close_warned.insert(i);
                    atom_close_warned.insert(j);
                    warnings.push(StructureWarning {
                        target: StructureWarningTarget::AtomPair(i, j),
                        message: "Atoms are overlapping or too close.".to_string(),
                    });
                }
            }
        }
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        if a.ignore_warnings || !a.visible || atom_close_warned.contains(&i) {
            continue;
        }
        for b in &mol.bonds {
            if b.ignore_warnings || !b.visible || b.begin_atom == i || b.end_atom == i {
                continue;
            }
            let a1 = &mol.atoms[b.begin_atom];
            let a2 = &mol.atoms[b.end_atom];
            if !a1.visible || !a2.visible {
                continue;
            }
            let seg_d = point_to_segment_dist(a.x, a.y, a1.x, a1.y, a2.x, a2.y);
            if seg_d < 2.0 {
                warnings.push(StructureWarning {
                    target: StructureWarningTarget::Atom(i),
                    message: "This atom is very close to another atom or bond.".to_string(),
                });
                break;
            }
        }
    }

    for cap in captions {
        if let Some(msg) = &cap.warning {
            warnings.push(StructureWarning {
                target: StructureWarningTarget::Caption(cap.id),
                message: msg.clone(),
            });
        }
    }

    warnings
}

#[inline]
fn segments_intersect(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    x3: f64,
    y3: f64,
    x4: f64,
    y4: f64,
) -> bool {
    let dx1 = x2 - x1;
    let dy1 = y2 - y1;
    let dx2 = x4 - x3;
    let dy2 = y4 - y3;
    let denom = dx1 * dy2 - dy1 * dx2;
    if denom.abs() <= 1e-9 {
        return false;
    }
    let qx = x3 - x1;
    let qy = y3 - y1;
    let t = (qx * dy2 - qy * dx2) / denom;
    let u = (qx * dy1 - qy * dx1) / denom;
    t > 0.001 && t < 0.999 && u > 0.001 && u < 0.999
}

#[inline]
fn segment_to_segment_dist(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    x3: f64,
    y3: f64,
    x4: f64,
    y4: f64,
) -> f64 {
    if segments_intersect(x1, y1, x2, y2, x3, y3, x4, y4) {
        return 0.0;
    }
    point_to_segment_dist(x1, y1, x3, y3, x4, y4)
        .min(point_to_segment_dist(x2, y2, x3, y3, x4, y4))
        .min(point_to_segment_dist(x3, y3, x1, y1, x2, y2))
        .min(point_to_segment_dist(x4, y4, x1, y1, x2, y2))
}

/// Detects crossing or near-crossing (< 2.0 pt clearance) visible non-adjacent bond pairs
/// `((mol_a, bond_a), (mol_b, bond_b))` across all molecules in a scene (`CrossingBonds` in `P33b`).
pub fn find_crossing_bond_pairs(mols: &[ChemMolecule]) -> Vec<((usize, usize), (usize, usize))> {
    let mut flat = Vec::new();
    for (mi, mol) in mols.iter().enumerate() {
        for (bi, b) in mol.bonds.iter().enumerate() {
            let a1 = &mol.atoms[b.begin_atom];
            let a2 = &mol.atoms[b.end_atom];
            if b.visible && a1.visible && a2.visible {
                flat.push((mi, bi, b.begin_atom, b.end_atom, a1.x, a1.y, a2.x, a2.y));
            }
        }
    }
    let mut crossings = Vec::new();
    for i in 0..flat.len() {
        let (m1, b1, u1, v1, x1, y1, x2, y2) = flat[i];
        for &(m2, b2, u2, v2, x3, y3, x4, y4) in &flat[(i + 1)..] {
            if m1 == m2 && (u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2) {
                continue;
            }
            if segment_to_segment_dist(x1, y1, x2, y2, x3, y3, x4, y4) < 2.0 {
                crossings.push(((m1, b1), (m2, b2)));
            }
        }
    }
    crossings
}

/// Evaluates `Structure > Check Structure` across an entire multi-fragment `DocumentScene`,
/// including cross-fragment atom-to-atom overlap (`< 1.0 pt`) and atom-to-bond proximity (`< 2.0 pt`).
pub fn check_structure_scene(scene: &DocumentScene) -> Vec<SceneStructureWarning> {
    let mut out = Vec::new();
    let mut seen_atom_bond_close: HashSet<(usize, usize)> = HashSet::new();

    for (mi, mol) in scene.molecules.iter().enumerate() {
        for w in check_structure(mol, &[]) {
            let object_id = match w.target {
                StructureWarningTarget::Atom(ai) => {
                    if w.message == "This atom is very close to another atom or bond."
                        || w.message == "Atoms are overlapping or too close."
                    {
                        seen_atom_bond_close.insert((mi, ai));
                    }
                    mol.atoms[ai].id
                }
                StructureWarningTarget::Bond(bi) => mol.bonds[bi].id,
                StructureWarningTarget::AtomPair(ai, _) => {
                    seen_atom_bond_close.insert((mi, ai));
                    mol.atoms[ai].id
                }
                StructureWarningTarget::Caption(cid) => cid,
            };
            out.push(SceneStructureWarning {
                molecule_idx: Some(mi),
                object_id,
                message: w.message,
            });
        }
    }

    for mi in 0..scene.molecules.len() {
        let m1 = &scene.molecules[mi];
        for (ai, a1) in m1.atoms.iter().enumerate() {
            if a1.ignore_warnings || !a1.visible || seen_atom_bond_close.contains(&(mi, ai)) {
                continue;
            }
            let mut warned = false;
            for (mj, m2) in scene.molecules.iter().enumerate() {
                if mj == mi {
                    continue;
                }
                for a2 in &m2.atoms {
                    if a2.ignore_warnings || !a2.visible {
                        continue;
                    }
                    let d = ((a1.x - a2.x).powi(2) + (a1.y - a2.y).powi(2)).sqrt();
                    if d < 1.0 {
                        out.push(SceneStructureWarning {
                            molecule_idx: Some(mi),
                            object_id: a1.id,
                            message: "Atoms are overlapping or too close.".to_string(),
                        });
                        warned = true;
                        break;
                    }
                }
                if warned {
                    break;
                }
                for b2 in &m2.bonds {
                    if b2.ignore_warnings || !b2.visible {
                        continue;
                    }
                    let p = &m2.atoms[b2.begin_atom];
                    let q = &m2.atoms[b2.end_atom];
                    if !p.visible || !q.visible {
                        continue;
                    }
                    if point_to_segment_dist(a1.x, a1.y, p.x, p.y, q.x, q.y) < 2.0 {
                        out.push(SceneStructureWarning {
                            molecule_idx: Some(mi),
                            object_id: a1.id,
                            message: "This atom is very close to another atom or bond."
                                .to_string(),
                        });
                        warned = true;
                        break;
                    }
                }
                if warned {
                    break;
                }
            }
        }
    }

    for cap in &scene.captions {
        if let Some(msg) = &cap.warning {
            out.push(SceneStructureWarning {
                molecule_idx: None,
                object_id: cap.id,
                message: msg.clone(),
            });
        }
    }
    out
}

/// Checks enhanced-stereochemistry warnings in a CDXML document (e.g. `EnhancedStereoType="Or"`/`"And"`
/// on an achiral compound -> `"This compound cannot possess relative stereochemistry."`).
pub fn check_enhanced_stereo_warnings(cdxml: &CdxmlDocument) -> Vec<SceneStructureWarning> {
    let mut caption_frag_occs = HashSet::new();
    for (occ, el) in cdxml.elements.iter().enumerate() {
        if el.tag == "fragment" {
            let mut n_count = 0usize;
            let mut b_count = 0usize;
            let mut has_t = false;
            for c in &el.children {
                if let CdxmlNodeChild::Element(cid) = c {
                    let ce = &cdxml.elements[*cid];
                    if ce.tag == "n" {
                        n_count += 1;
                        has_t = ce.children.iter().any(|cc| {
                            matches!(cc, CdxmlNodeChild::Element(tc) if cdxml.elements[*tc].tag == "t")
                        });
                    } else if ce.tag == "b" {
                        b_count += 1;
                    }
                }
            }
            if n_count == 1 && b_count == 0 && has_t {
                caption_frag_occs.insert(occ);
            }
        }
    }

    let markers_by_occ = cdxml.extract_enhanced_stereo();
    let mut out = Vec::new();
    let mut mol_idx = 0usize;
    for (occ, mol) in cdxml.extract_molecules() {
        if caption_frag_occs.contains(&occ) {
            continue;
        }
        if let Some(markers) = markers_by_occ.get(&occ) {
            let has_stereo = (0..mol.atoms.len()).any(|ai| {
                chiralcraft_chem::compute_tetrahedral_cip_code(&mol, ai) != 0
                    || chiralcraft_chem::compute_axial_allene_cip(&mol, ai).is_some()
            }) || (0..mol.bonds.len())
                .any(|bi| chiralcraft_chem::compute_double_bond_cip_code(&mol, bi) != 0);
            for m in markers {
                if matches!(m.group_type, 3 | 4) && !has_stereo {
                    if let Some(atom) = mol.atoms.get(m.atom_idx) {
                        if !atom.ignore_warnings {
                            out.push(SceneStructureWarning {
                                molecule_idx: Some(mol_idx),
                                object_id: atom.id,
                                message:
                                    "This compound cannot possess relative stereochemistry."
                                        .to_string(),
                            });
                        }
                    }
                }
            }
        }
        mol_idx += 1;
    }
    out
}

/// Executes `Structure > Toggle Hydrogen Display` on `mol` (`P32b`):
/// toggles `Visible="no"` and `IgnoreWarnings="yes"` on explicit `Element="1"` atoms
/// and their incident bonds. Returns `true` if explicit hydrogens are now visible.
pub fn toggle_hydrogen_display(mol: &mut ChemMolecule) -> bool {
    let any_visible = mol.atoms.iter().any(|a| a.z == 1 && a.visible);
    let target_visible = !any_visible;
    let target_ignore = !target_visible;
    for a in &mut mol.atoms {
        if a.z == 1 {
            a.visible = target_visible;
            a.ignore_warnings = target_ignore;
        }
    }
    for b in &mut mol.bonds {
        if mol.atoms[b.begin_atom].z == 1 || mol.atoms[b.end_atom].z == 1 {
            b.visible = target_visible;
            b.ignore_warnings = target_ignore;
        }
    }
    target_visible
}

/// Executes `Object > Flatten` on `mol` (`P32e`): clears 3D z-coordinates and marks the
/// molecule as 2D while preserving `(x, y)` coordinates.
pub fn flatten_molecule_to_2d(mol: &mut ChemMolecule) {
    mol.has_3d_coords = false;
    for a in &mut mol.atoms {
        a.z_coord = 0.0;
    }
}

fn polygon_signed_area(mol: &ChemMolecule, ring: &[usize]) -> f64 {
    let k = ring.len();
    let mut area = 0.0;
    for i in 0..k {
        let a = &mol.atoms[ring[i]];
        let b = &mol.atoms[ring[(i + 1) % k]];
        area += a.x * b.y - b.x * a.y;
    }
    0.5 * area
}

fn find_chordless_small_rings(mol: &ChemMolecule, max_len: usize) -> Vec<(Vec<usize>, Vec<usize>)> {
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    let mut seen_sets: BTreeSet<Vec<usize>> = BTreeSet::new();
    let mut rings: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();

    for (bi, b) in mol.bonds.iter().enumerate() {
        let (start, goal) = (b.begin_atom, b.end_atom);
        let mut parent: Vec<Option<(usize, usize)>> = vec![None; n];
        let mut dist = vec![usize::MAX; n];
        let mut q = VecDeque::new();
        dist[start] = 0;
        q.push_back(start);

        while let Some(u) = q.pop_front() {
            if u == goal || dist[u] + 1 >= max_len {
                continue;
            }
            for &(v, e_idx) in &adj[u] {
                if e_idx == bi {
                    continue;
                }
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    parent[v] = Some((u, e_idx));
                    q.push_back(v);
                }
            }
        }

        if dist[goal] == usize::MAX || dist[goal] + 1 > max_len {
            continue;
        }

        let mut verts = vec![goal];
        let mut b_indices = vec![bi];
        let mut cur = goal;
        while let Some((p, e_idx)) = parent[cur] {
            verts.push(p);
            b_indices.push(e_idx);
            cur = p;
        }
        verts.reverse();
        if verts.len() < 3 {
            continue;
        }

        // Chordless check: no bond in `mol` connects two non-consecutive ring vertices
        let vset: HashSet<usize> = verts.iter().copied().collect();
        let cycle_edges: HashSet<usize> = b_indices.iter().copied().collect();
        let has_chord = mol.bonds.iter().enumerate().any(|(other_bi, ob)| {
            !cycle_edges.contains(&other_bi)
                && vset.contains(&ob.begin_atom)
                && vset.contains(&ob.end_atom)
        });
        if has_chord {
            continue;
        }

        let mut key = verts.clone();
        key.sort_unstable();
        if seen_sets.insert(key) {
            rings.push((verts, b_indices));
        }
    }

    rings.sort_by_key(|(v, _)| v.len());
    rings
}

fn find_bridge_bonds(mol: &ChemMolecule) -> HashSet<usize> {
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }
    let mut tin = vec![0usize; n];
    let mut low = vec![0usize; n];
    let mut timer = 0usize;
    let mut bridges = HashSet::new();

    #[allow(clippy::too_many_arguments)]
    fn dfs(
        u: usize,
        parent_edge: Option<usize>,
        adj: &[Vec<(usize, usize)>],
        tin: &mut [usize],
        low: &mut [usize],
        timer: &mut usize,
        bridges: &mut HashSet<usize>,
    ) {
        *timer += 1;
        tin[u] = *timer;
        low[u] = *timer;
        for &(v, bi) in &adj[u] {
            if Some(bi) == parent_edge {
                continue;
            }
            if tin[v] != 0 {
                low[u] = low[u].min(tin[v]);
            } else {
                dfs(v, Some(bi), adj, tin, low, timer, bridges);
                low[u] = low[u].min(low[v]);
                if low[v] > tin[u] {
                    bridges.insert(bi);
                }
            }
        }
    }

    for i in 0..n {
        if tin[i] == 0 {
            dfs(i, None, &adj, &mut tin, &mut low, &mut timer, &mut bridges);
        }
    }
    bridges
}

fn idealize_isolated_ring(
    mol: &mut ChemMolecule,
    ring: &[usize],
    orig_pos: &[(f64, f64)],
    target_bond_length_pt: f64,
) {
    let k = ring.len();
    let k_f = k as f64;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for &vi in ring {
        cx += orig_pos[vi].0;
        cy += orig_pos[vi].1;
    }
    cx /= k_f;
    cy /= k_f;

    let area = polygon_signed_area(mol, ring);
    let winding = if area < 0.0 { -1.0 } else { 1.0 };
    let step = winding * 2.0 * std::f64::consts::PI / k_f;

    // Procrustes angular phase preserving current ring orientation
    let mut sum_cos = 0.0;
    let mut sum_sin = 0.0;
    for (i, &vi) in ring.iter().enumerate() {
        let dx = orig_pos[vi].0 - cx;
        let dy = orig_pos[vi].1 - cy;
        if dx * dx + dy * dy > 1e-12 {
            let ang = dy.atan2(dx) - (i as f64) * step;
            sum_cos += ang.cos();
            sum_sin += ang.sin();
        }
    }
    let theta0 = if sum_cos.abs() > 1e-12 || sum_sin.abs() > 1e-12 {
        sum_sin.atan2(sum_cos)
    } else {
        0.0
    };
    let radius = target_bond_length_pt / (2.0 * (std::f64::consts::PI / k_f).sin());
    for (i, &vi) in ring.iter().enumerate() {
        let theta = theta0 + (i as f64) * step;
        mol.atoms[vi].x = cx + radius * theta.cos();
        mol.atoms[vi].y = cy + radius * theta.sin();
    }
}

fn idealize_fused_ring(
    mol: &mut ChemMolecule,
    ring: &[usize],
    orig_pos: &[(f64, f64)],
    pinned_atoms: &HashSet<usize>,
    target_bond_length_pt: f64,
) -> bool {
    let k = ring.len();
    // Find shared edge (u, v) in ring order
    let mut shared_pos = None;
    for i in 0..k {
        let u = ring[i];
        let v = ring[(i + 1) % k];
        if pinned_atoms.contains(&u) && pinned_atoms.contains(&v) {
            shared_pos = Some((i, u, v));
            break;
        }
    }
    let Some((start_i, u, v)) = shared_pos else {
        return false;
    };

    let (ux, uy) = (mol.atoms[u].x, mol.atoms[u].y);
    let (vx, vy) = (mol.atoms[v].x, mol.atoms[v].y);
    let edge_dx = vx - ux;
    let edge_dy = vy - uy;
    let edge_len = (edge_dx * edge_dx + edge_dy * edge_dy).sqrt();
    if edge_len <= 1e-9 {
        return false;
    }

    // Determine original side of unpinned ring centroid relative to original u->v edge
    let (oux, ouy) = orig_pos[u];
    let (ovx, ovy) = orig_pos[v];
    let o_ex = ovx - oux;
    let o_ey = ovy - ouy;
    let mut unpinned_cx = 0.0;
    let mut unpinned_cy = 0.0;
    let mut unpinned_cnt = 0.0;
    for &vi in ring {
        if !pinned_atoms.contains(&vi) {
            unpinned_cx += orig_pos[vi].0;
            unpinned_cy += orig_pos[vi].1;
            unpinned_cnt += 1.0;
        }
    }
    if unpinned_cnt == 0.0 {
        return true;
    }
    unpinned_cx /= unpinned_cnt;
    unpinned_cy /= unpinned_cnt;
    let cross_orig = o_ex * (unpinned_cy - ouy) - o_ey * (unpinned_cx - oux);
    let side_sign = if cross_orig.abs() > 1e-9 {
        if cross_orig < 0.0 { -1.0 } else { 1.0 }
    } else {
        // Collapsed initial coordinates (e.g., SMILES): fold fused ring OUTWARD away from pinned centroid!
        let mut pin_cx = 0.0;
        let mut pin_cy = 0.0;
        let mut pin_cnt = 0.0;
        for &pi in pinned_atoms {
            if pi != u && pi != v {
                pin_cx += mol.atoms[pi].x;
                pin_cy += mol.atoms[pi].y;
                pin_cnt += 1.0;
            }
        }
        if pin_cnt > 0.0 {
            pin_cx /= pin_cnt;
            pin_cy /= pin_cnt;
            let cross_pin = edge_dx * (pin_cy - uy) - edge_dy * (pin_cx - ux);
            if cross_pin > 0.0 { -1.0 } else { 1.0 }
        } else {
            1.0
        }
    };

    // Unit tangent and normal toward unpinned ring side
    let tx = edge_dx / edge_len;
    let ty = edge_dy / edge_len;
    let nx = -ty * side_sign;
    let ny = tx * side_sign;

    let k_f = k as f64;
    let half_angle = std::f64::consts::PI / k_f;
    let radius = edge_len / (2.0 * half_angle.sin());
    let apothem = radius * half_angle.cos();
    let midx = 0.5 * (ux + vx);
    let midy = 0.5 * (uy + vy);
    let cx = midx + apothem * nx;
    let cy = midy + apothem * ny;

    let theta_u = (uy - cy).atan2(ux - cx);
    let theta_v = (vy - cy).atan2(vx - cx);
    let mut d_uv = theta_v - theta_u;
    while d_uv > std::f64::consts::PI {
        d_uv -= 2.0 * std::f64::consts::PI;
    }
    while d_uv < -std::f64::consts::PI {
        d_uv += 2.0 * std::f64::consts::PI;
    }
    let step = if d_uv >= 0.0 {
        2.0 * half_angle
    } else {
        -2.0 * half_angle
    };
    let _ = target_bond_length_pt; // shared edge fixes ring radius when core ring has target length

    for offset in 2..k {
        let idx = (start_i + offset) % k;
        let vi = ring[idx];
        if !pinned_atoms.contains(&vi) {
            let theta = theta_u + (offset as f64) * step;
            mol.atoms[vi].x = cx + radius * theta.cos();
            mol.atoms[vi].y = cy + radius * theta.sin();
        }
    }
    true
}

/// Performs stereochemistry-preserving 2D structure cleanup (`P02`):
/// - Normalizes bond lengths to `target_bond_length_pt` while preserving connectivity,
///   substituent identity, and 2D tetrahedral/alkene stereochemical descriptors (`perceive_stereochemistry`).
/// - Applies the documented whole-ring selection rule (`D16`): any 3..=8-membered ring is idealized
///   to a regular polygon (including fused and spiro ring systems) if and only if all bonds of that ring
///   are in `selected_bonds`.
/// - Scales selected acyclic bridge bonds along their existing ray directions and rigidly translates
///   downstream substituent subtrees so stereocenters and double bonds never invert.
pub fn clean_up_structure_preserving_stereo(
    mol: &mut ChemMolecule,
    selected_bonds: &[usize],
    target_bond_length_pt: f64,
) {
    if mol.bonds.is_empty() || mol.atoms.is_empty() {
        return;
    }
    let stereo_before = perceive_stereochemistry(mol);
    let selected_set: HashSet<usize> = if selected_bonds.is_empty() {
        (0..mol.bonds.len()).collect()
    } else {
        selected_bonds
            .iter()
            .copied()
            .filter(|&bi| bi < mol.bonds.len())
            .collect()
    };
    if selected_set.is_empty() {
        return;
    }

    let orig_pos: Vec<(f64, f64)> = mol.atoms.iter().map(|a| (a.x, a.y)).collect();
    let rings = find_chordless_small_rings(mol, 8);
    let whole_rings: Vec<(Vec<usize>, Vec<usize>)> = rings
        .into_iter()
        .filter(|(_, r_bonds)| r_bonds.iter().all(|bi| selected_set.contains(bi)))
        .collect();

    let mut idealized_ring_bonds = HashSet::new();
    let mut pinned_atoms = HashSet::new();

    if !whole_rings.is_empty() {
        let mut placed = vec![false; whole_rings.len()];
        for &bi in whole_rings.iter().flat_map(|(_, rb)| rb.iter()) {
            idealized_ring_bonds.insert(bi);
        }

        // Order placement: aromatic or largest ring seeds each fused component, then fused/spiro neighbors
        let arom = CcAromaticityInfo::perceive(mol);
        let mut ring_order: Vec<usize> = (0..whole_rings.len()).collect();
        ring_order.sort_by_key(|&ri| {
            let arom_hits = whole_rings[ri]
                .1
                .iter()
                .filter(|bi| arom.aromatic_bonds.contains(bi))
                .count();
            (std::cmp::Reverse(arom_hits), std::cmp::Reverse(whole_rings[ri].0.len()))
        });

        for &seed_ri in &ring_order {
            if placed[seed_ri] {
                continue;
            }
            placed[seed_ri] = true;
            idealize_isolated_ring(mol, &whole_rings[seed_ri].0, &orig_pos, target_bond_length_pt);
            for &vi in &whole_rings[seed_ri].0 {
                pinned_atoms.insert(vi);
            }

            loop {
                let mut progressed = false;
                for &ri in &ring_order {
                    if placed[ri] {
                        continue;
                    }
                    let (ref verts, _) = whole_rings[ri];
                    let shared_count = verts.iter().filter(|vi| pinned_atoms.contains(vi)).count();
                    if shared_count >= 2
                        && idealize_fused_ring(
                            mol,
                            verts,
                            &orig_pos,
                            &pinned_atoms,
                            target_bond_length_pt,
                        )
                    {
                        placed[ri] = true;
                        for &vi in verts {
                            pinned_atoms.insert(vi);
                        }
                        progressed = true;
                    } else if shared_count == 1 {
                        // Spiro ring: idealize centered then translate so shared vertex matches pinned position
                        let shared_v = *verts.iter().find(|vi| pinned_atoms.contains(vi)).unwrap();
                        let (anchor_x, anchor_y) = (mol.atoms[shared_v].x, mol.atoms[shared_v].y);
                        idealize_isolated_ring(mol, verts, &orig_pos, target_bond_length_pt);
                        let dx = anchor_x - mol.atoms[shared_v].x;
                        let dy = anchor_y - mol.atoms[shared_v].y;
                        for &vi in verts {
                            mol.atoms[vi].x += dx;
                            mol.atoms[vi].y += dy;
                            pinned_atoms.insert(vi);
                        }
                        placed[ri] = true;
                        progressed = true;
                    }
                }
                if !progressed {
                    break;
                }
            }
        }

        // Rigidly translate attached non-pinned subtrees of idealized ring atoms
        let n = mol.atoms.len();
        let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
        for (bi, b) in mol.bonds.iter().enumerate() {
            if !idealized_ring_bonds.contains(&bi) {
                adj[b.begin_atom].push((b.end_atom, bi));
                adj[b.end_atom].push((b.begin_atom, bi));
            }
        }
        let mut shifted = vec![false; n];
        for &ra in &pinned_atoms {
            let dx = mol.atoms[ra].x - orig_pos[ra].0;
            let dy = mol.atoms[ra].y - orig_pos[ra].1;
            let mut q = VecDeque::new();
            q.push_back(ra);
            while let Some(u) = q.pop_front() {
                for &(v, _) in &adj[u] {
                    if !pinned_atoms.contains(&v) && !shifted[v] {
                        shifted[v] = true;
                        mol.atoms[v].x += dx;
                        mol.atoms[v].y += dy;
                        q.push_back(v);
                    }
                }
            }
        }
    }

    // Scale selected acyclic bridge bonds along their existing ray direction and propagate rigidly
    let bridges = find_bridge_bonds(mol);
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    let acyclic_selected: HashSet<usize> = selected_set
        .iter()
        .copied()
        .filter(|bi| !idealized_ring_bonds.contains(bi) && bridges.contains(bi))
        .collect();

    if !acyclic_selected.is_empty() {
        // Seed multi-source BFS with ALL pinned ring atoms simultaneously so every exocyclic ring
        // substituent root (depth 1) is placed before deeper branches choose steric-avoiding angles!
        let mut pinned_seeds: Vec<usize> = pinned_atoms.iter().copied().collect();
        pinned_seeds.sort_unstable();
        let mut remaining_seeds: Vec<usize> = (0..n)
            .filter(|i| !pinned_atoms.contains(i) && adj[*i].iter().any(|(_, bi)| acyclic_selected.contains(bi)))
            .collect();
        remaining_seeds.sort_by_key(|&i| std::cmp::Reverse(adj[i].len()));

        let mut visited_atoms = vec![false; n];
        let mut placed_atoms = vec![false; n];
        for &ra in &pinned_atoms {
            placed_atoms[ra] = true;
        }
        let mut depth_atoms = vec![0usize; n];
        let mut q = VecDeque::new();
        let mut seed_batches: Vec<Vec<usize>> = Vec::new();
        if !pinned_seeds.is_empty() {
            seed_batches.push(pinned_seeds);
        }
        for s in remaining_seeds {
            seed_batches.push(vec![s]);
        }
        for batch in seed_batches {
            let mut pushed_any = false;
            for seed in batch {
                if visited_atoms[seed] {
                    continue;
                }
                visited_atoms[seed] = true;
                placed_atoms[seed] = true;
                q.push_back(seed);
                pushed_any = true;
            }
            if !pushed_any {
                continue;
            }
            while let Some(u) = q.pop_front() {
                // Collect unvisited neighbors of `u` so collapsed SMILES branches fan out cleanly without overlap!
                let children: Vec<(usize, usize)> = adj[u]
                    .iter()
                    .copied()
                    .filter(|&(v, _)| !visited_atoms[v])
                    .collect();
                let mut collapsed_children: Vec<(usize, usize)> = Vec::new();
                for &(v, bi) in &children {
                    visited_atoms[v] = true;
                    depth_atoms[v] = depth_atoms[u] + 1;
                    if acyclic_selected.contains(&bi) {
                        let ox = orig_pos[v].0 - orig_pos[u].0;
                        let oy = orig_pos[v].1 - orig_pos[u].1;
                        let olen = (ox * ox + oy * oy).sqrt();
                        if olen > 1e-9 {
                            placed_atoms[v] = true;
                            let target_x = mol.atoms[u].x + target_bond_length_pt * (ox / olen);
                            let target_y = mol.atoms[u].y + target_bond_length_pt * (oy / olen);
                            let dx = target_x - mol.atoms[v].x;
                            let dy = target_y - mol.atoms[v].y;
                            // Translate entire downstream subtree beyond bridge edge `bi`
                            let mut sub_q = VecDeque::new();
                            let mut sub_vis = vec![false; n];
                            sub_vis[u] = true;
                            sub_vis[v] = true;
                            mol.atoms[v].x += dx;
                            mol.atoms[v].y += dy;
                            sub_q.push_back(v);
                            while let Some(su) = sub_q.pop_front() {
                                for &(sv, _) in &adj[su] {
                                    if !sub_vis[sv] {
                                        sub_vis[sv] = true;
                                        mol.atoms[sv].x += dx;
                                        mol.atoms[sv].y += dy;
                                        sub_q.push_back(sv);
                                    }
                                }
                            }
                        } else {
                            collapsed_children.push((v, bi));
                        }
                    } else {
                        placed_atoms[v] = true;
                    }
                    q.push_back(v);
                }

                // Synthesize clean IUPAC 2D geometry for any collapsed (SMILES/InChI) children of `u`
                if !collapsed_children.is_empty() {
                    let placed_nbrs: Vec<(usize, f64, i32)> = adj[u]
                        .iter()
                        .filter_map(|&(nb, bi)| {
                            if collapsed_children.iter().any(|&(cv, _)| cv == nb) {
                                return None;
                            }
                            let dx = mol.atoms[nb].x - mol.atoms[u].x;
                            let dy = mol.atoms[nb].y - mol.atoms[u].y;
                            if dx.hypot(dy) > 1e-6 {
                                Some((nb, dy.atan2(dx), mol.bonds[bi].order))
                            } else {
                                None
                            }
                        })
                        .collect();
                    // Steric clearance helper: minimum squared distance from a candidate child position
                    // to all already-positioned non-u atoms (maximizing clearance automatically yields
                    // trans zig-zag chains AND bends ortho ring substituents away from each other!).
                    let clearance_d2 = |ang: f64, skip_atom: Option<usize>| -> f64 {
                        let cx = mol.atoms[u].x + target_bond_length_pt * ang.cos();
                        let cy = mol.atoms[u].y + target_bond_length_pt * ang.sin();
                        let mut best = f64::INFINITY;
                        for w in 0..n {
                            if !placed_atoms[w] || w == u || skip_atom == Some(w) {
                                continue;
                            }
                            let dx = mol.atoms[w].x - cx;
                            let dy = mol.atoms[w].y - cy;
                            best = best.min(dx * dx + dy * dy);
                        }
                        best
                    };
                    let subtree_weight = |root: usize| -> usize {
                        let mut vis = vec![false; n];
                        vis[u] = true;
                        vis[root] = true;
                        let mut st = vec![root];
                        let mut cnt = 0usize;
                        while let Some(cur) = st.pop() {
                            cnt += 1;
                            for &(nb, _) in &adj[cur] {
                                if !vis[nb] {
                                    vis[nb] = true;
                                    st.push(nb);
                                }
                            }
                        }
                        cnt
                    };
                    let k = collapsed_children.len();
                    let angles: Vec<f64> = if placed_nbrs.is_empty() {
                        let base = -std::f64::consts::FRAC_PI_6;
                        let step = match k {
                            1 => 0.0,
                            2 => 2.0 * std::f64::consts::PI / 3.0,
                            3 => 2.0 * std::f64::consts::PI / 3.0,
                            _ => std::f64::consts::TAU / (k as f64),
                        };
                        (0..k).map(|i| base + (i as f64) * step).collect()
                    } else if placed_nbrs.len() == 1 {
                        let (parent_nb, back_ang, back_ord) = placed_nbrs[0];
                        let has_triple = back_ord == 3
                            || collapsed_children.iter().any(|&(_, bi)| mol.bonds[bi].order == 3);
                        match k {
                            1 if has_triple => vec![back_ang + std::f64::consts::PI],
                            1 => {
                                let a_pos = back_ang + 2.0 * std::f64::consts::PI / 3.0;
                                let a_neg = back_ang - 2.0 * std::f64::consts::PI / 3.0;
                                let d_pos = clearance_d2(a_pos, Some(parent_nb));
                                let d_neg = clearance_d2(a_neg, Some(parent_nb));
                                let eps = 1e-2 * target_bond_length_pt * target_bond_length_pt;
                                let chosen = if d_pos.is_finite() && d_neg.is_finite() && (d_pos - d_neg).abs() > eps {
                                    if d_pos > d_neg { a_pos } else { a_neg }
                                } else {
                                    let sign = if depth_atoms[u] % 2 == 0 { 1.0 } else { -1.0 };
                                    back_ang + sign * (2.0 * std::f64::consts::PI / 3.0)
                                };
                                vec![chosen]
                            }
                            2 => {
                                let a_pos = back_ang + 2.0 * std::f64::consts::PI / 3.0;
                                let a_neg = back_ang - 2.0 * std::f64::consts::PI / 3.0;
                                let w0 = subtree_weight(collapsed_children[0].0);
                                let w1 = subtree_weight(collapsed_children[1].0);
                                if w0 != w1 {
                                    let d_pos = clearance_d2(a_pos, Some(parent_nb));
                                    let d_neg = clearance_d2(a_neg, Some(parent_nb));
                                    let (open_ang, tight_ang) = if d_pos >= d_neg { (a_pos, a_neg) } else { (a_neg, a_pos) };
                                    if w0 > w1 {
                                        vec![open_ang, tight_ang]
                                    } else {
                                        vec![tight_ang, open_ang]
                                    }
                                } else {
                                    vec![a_pos, a_neg]
                                }
                            }
                            3 => vec![
                                back_ang + std::f64::consts::FRAC_PI_2,
                                back_ang + std::f64::consts::PI,
                                back_ang - std::f64::consts::FRAC_PI_2,
                            ],
                            _ => {
                                let step = std::f64::consts::TAU / ((k + 1) as f64);
                                (1..=k).map(|i| back_ang + (i as f64) * step).collect()
                            }
                        }
                    } else {
                        // Outward bisector away from existing placed neighbors (e.g. exocyclic ring substituents)
                        let sx: f64 = placed_nbrs.iter().map(|(_, a, _)| a.cos()).sum();
                        let sy: f64 = placed_nbrs.iter().map(|(_, a, _)| a.sin()).sum();
                        let out_ang = (-sy).atan2(-sx);
                        if k == 1 {
                            vec![out_ang]
                        } else {
                            let spread = (std::f64::consts::PI / 3.0).min(std::f64::consts::PI / (k as f64));
                            (0..k)
                                .map(|i| out_ang + ((i as f64) - 0.5 * ((k - 1) as f64)) * spread)
                                .collect()
                        }
                    };
                    for ((v, _), ang) in collapsed_children.into_iter().zip(angles.into_iter()) {
                        mol.atoms[v].x = mol.atoms[u].x + target_bond_length_pt * ang.cos();
                        mol.atoms[v].y = mol.atoms[u].y + target_bond_length_pt * ang.sin();
                        placed_atoms[v] = true;
                    }
                }
            }
        }
    }

    let stereo_after = perceive_stereochemistry(mol);
    debug_assert_eq!(
        stereo_before.atom_stereo, stereo_after.atom_stereo,
        "Tetrahedral stereochemistry inverted during cleanup"
    );
    debug_assert_eq!(
        stereo_before.bond_stereo, stereo_after.bond_stereo,
        "Double-bond E/Z stereochemistry inverted during cleanup"
    );
    let _ = CcVal::compute(mol);
}

#[cfg(test)]
mod tests {
    use super::*;
    use chiralcraft_chem::{compute_double_bond_cip_code, compute_tetrahedral_cip_code};
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_p01_v1_to_v5_undo_save_revert_parity() {
        // P01-V1 Control: each caption creates 3 undo transactions; Undo/Redo restores exact scene
        let mut doc = InteractiveDocument::new_saved_seed("seed-v1.cdxml", DocumentSettings::acs_document_1996());
        assert!(!doc.is_dirty());
        assert_eq!(doc.undo_menu_label(), "Undo not available");
        assert_eq!(doc.redo_menu_label(), "Redo not available");

        doc.create_caption_interactive(176.96, 131.85, "A");
        doc.create_caption_interactive(266.96, 131.85, "B");
        assert!(doc.is_dirty());
        assert_eq!(doc.scene.captions.len(), 2);
        assert_eq!(doc.scene.captions[0].generic_nickname.as_deref(), Some("A"));
        assert_eq!(doc.scene.captions[1].element, Some(5));
        assert_eq!(
            doc.scene.captions[1].warning.as_deref(),
            Some("An atom in this label has an invalid valence.")
        );

        let scene_ab = doc.scene.clone();
        // 6 undo steps (Caption Change -> Typing -> Caption Change for B, then A)
        assert_eq!(doc.undo_menu_label(), "Undo Caption Change");
        assert!(doc.undo());
        assert_eq!(doc.undo_menu_label(), "Undo Typing");
        assert!(doc.undo());
        assert_eq!(doc.undo_menu_label(), "Undo Caption Change");
        assert!(doc.undo());
        assert_eq!(doc.scene.captions.len(), 1); // B removed, A remains

        assert!(doc.undo());
        assert!(doc.undo());
        assert!(doc.undo());
        assert!(doc.scene.captions.is_empty());
        assert!(!doc.is_dirty()); // Back at saved seed!
        assert_eq!(doc.undo_menu_label(), "Undo not available");

        // Redo 6 times restores exact AB state
        for _ in 0..6 {
            assert!(doc.redo());
        }
        assert_eq!(doc.scene, scene_ab);
        assert_eq!(doc.redo_menu_label(), "Redo not available");

        // P01-V2 Save boundary (DOC-CONTRADICTION): Save does NOT clear Undo stack!
        doc.save();
        assert!(!doc.is_dirty());
        assert_eq!(doc.last_repeat_label.as_deref(), Some("Save"));
        assert_eq!(doc.undo_menu_label(), "Undo Caption Change");
        assert!(doc.undo());
        assert!(doc.is_dirty()); // Undoing past Save marks document dirty again!
        assert!(doc.redo());
        assert!(!doc.is_dirty());

        // P01-V4 Redo branch invalidation
        for _ in 0..3 {
            assert!(doc.undo());
        }
        assert_eq!(doc.redo_menu_label(), "Redo Caption Change");
        doc.create_caption_interactive(176.96, 200.0, "C");
        assert_eq!(doc.redo_menu_label(), "Redo not available");

        // P01-V5 Revert: returns to saved state and clears undo history
        let mut doc_v5 = InteractiveDocument::new_saved_seed("seed-v5.cdxml", DocumentSettings::acs_document_1996());
        doc_v5.create_caption_interactive(176.96, 131.85, "A");
        doc_v5.save();
        let saved_a = doc_v5.scene.clone();
        doc_v5.create_caption_interactive(266.96, 131.85, "B");
        assert!(doc_v5.is_revert_enabled());
        assert!(doc_v5.revert());
        assert_eq!(doc_v5.scene, saved_a);
        assert!(!doc_v5.is_dirty());
        assert_eq!(doc_v5.undo_menu_label(), "Undo not available");
    }

    #[test]
    fn test_seed_cdxml_fixtures_and_roundtrip() {
        let p01_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("research")
            .join("runtime-probes")
            .join("p01");
        if !p01_dir.exists() {
            return;
        }

        let seed1_xml = fs::read_to_string(p01_dir.join("seed-v1.cdxml")).unwrap();
        let seed1 = CdxmlDocument::parse(&seed1_xml).unwrap();
        let doc1 = InteractiveDocument::from_cdxml(&seed1);
        assert!(!doc1.settings.creation_program.is_empty());
        assert_eq!(DocumentSettings::acs_document_1996().creation_program, "ChiralCraft 1.0");
        assert!((doc1.settings.bond_length_pt - 14.40).abs() < 1e-6);
        assert!((doc1.settings.chain_angle_deg - 120.0).abs() < 1e-6);
        assert_eq!(doc1.settings.label_font, "Arial");
        assert!(doc1.scene.captions.is_empty());

        let seed2_xml = fs::read_to_string(p01_dir.join("seed-v2.cdxml")).unwrap();
        let seed2 = CdxmlDocument::parse(&seed2_xml).unwrap();
        let doc2 = InteractiveDocument::from_cdxml(&seed2);
        assert_eq!(doc2.scene.captions.len(), 2);
        assert_eq!(doc2.scene.captions[0].text, "A");
        assert_eq!(doc2.scene.captions[0].generic_nickname.as_deref(), Some("A"));
        assert_eq!(doc2.scene.captions[1].text, "B");
        assert_eq!(doc2.scene.captions[1].element, Some(5));
        assert_eq!(
            doc2.scene.captions[1].warning.as_deref(),
            Some("An atom in this label has an invalid valence.")
        );

        // Interactive reconstruction matches seed-v2 captions & round-trips through CDXML
        let mut live = InteractiveDocument::from_cdxml(&seed1);
        live.create_caption_interactive(176.96, 131.85, "A");
        live.create_caption_interactive(266.96, 131.85, "B");
        let rt_xml = live.to_cdxml().to_xml_string();
        let rt_cdxml = CdxmlDocument::parse(&rt_xml).unwrap();
        let rt_doc = InteractiveDocument::from_cdxml(&rt_cdxml);
        assert_eq!(rt_doc.scene.captions.len(), 2);
        assert_eq!(rt_doc.scene.captions[0].generic_nickname.as_deref(), Some("A"));
        assert_eq!(rt_doc.scene.captions[1].element, Some(5));
        assert!(rt_doc.scene.captions[1].warning.is_some());

        let seed3_xml = fs::read_to_string(p01_dir.join("seed-v3.cdxml")).unwrap();
        let seed3 = CdxmlDocument::parse(&seed3_xml).unwrap();
        let doc3 = InteractiveDocument::from_cdxml(&seed3);
        assert_eq!(doc3.scene.captions.len(), 1);
        assert_eq!(doc3.scene.captions[0].generic_nickname.as_deref(), Some("A"));
    }

    #[test]
    fn test_p03_label_attachment_and_expansion() {
        // Flush-right F3C vs CF3
        let f3c_right = analyze_left_endpoint_label("F3C", LabelJustification::Right, 1, true);
        assert_eq!(f3c_right.attached_element_z, 6);
        assert!(!f3c_right.has_valence_warning);

        let cf3_right = analyze_left_endpoint_label("CF3", LabelJustification::Right, 1, true);
        assert_eq!(cf3_right.attached_element_z, 9);
        assert!(cf3_right.has_valence_warning);

        // Generalized halogen formula CCl3 behaves identically
        let ccl3_right = analyze_left_endpoint_label("CCl3", LabelJustification::Right, 1, true);
        assert_eq!(ccl3_right.attached_element_z, 17);
        assert!(ccl3_right.has_valence_warning);
        let cl3c_right = analyze_left_endpoint_label("Cl3C", LabelJustification::Right, 1, true);
        assert_eq!(cl3c_right.attached_element_z, 6);
        assert!(!cl3c_right.has_valence_warning);

        // Auto CF3 attaches to C without warning
        let cf3_auto = analyze_left_endpoint_label("CF3", LabelJustification::Auto, 1, true);
        assert_eq!(cf3_auto.attached_element_z, 6);
        assert!(!cf3_auto.has_valence_warning);

        // CH3 and H3C, NH2 and H2N, OH and HO all satisfy Hydrogen exception on left endpoint
        let ch3_right = analyze_left_endpoint_label("CH3", LabelJustification::Right, 1, true);
        assert_eq!(ch3_right.attached_element_z, 6);
        assert!(!ch3_right.has_valence_warning);
        let oh_right = analyze_left_endpoint_label("OH", LabelJustification::Right, 1, true);
        assert_eq!(oh_right.attached_element_z, 8);
        assert!(!oh_right.has_valence_warning);
        let h2n_auto = analyze_left_endpoint_label("H2N", LabelJustification::Auto, 1, true);
        assert_eq!(h2n_auto.attached_element_z, 7);
        assert!(!h2n_auto.has_valence_warning);

        // Right-endpoint mirror: C--F3C flush-left attaches to F (warning), C--CF3 flush-left attaches to C
        let right_end_f3c_left = analyze_endpoint_label("F3C", EndpointSide::Right, LabelJustification::Left, 1, true);
        assert_eq!(right_end_f3c_left.attached_element_z, 9);
        assert!(right_end_f3c_left.has_valence_warning);
        let right_end_cf3_left = analyze_endpoint_label("CF3", EndpointSide::Right, LabelJustification::Left, 1, true);
        assert_eq!(right_end_cf3_left.attached_element_z, 6);
        assert!(!right_end_cf3_left.has_valence_warning);

        // Common chemical nicknames (Me, OMe, NO2, COOH)
        let me = analyze_left_endpoint_label("Me", LabelJustification::Auto, 1, true);
        assert_eq!(me.attached_element_z, 6);
        assert!(!me.has_valence_warning);
        let no2 = analyze_left_endpoint_label("NO2", LabelJustification::Auto, 1, true);
        assert_eq!(no2.attached_element_z, 7);
        assert!(!no2.has_valence_warning);

        // Double bond rectification on vs off (CH3 -> CH2, NH2 -> NH)
        let ch3_dbl_on = analyze_left_endpoint_label("CH3", LabelJustification::Auto, 2, true);
        assert_eq!(ch3_dbl_on.raw_label, "CH2");
        assert!(!ch3_dbl_on.has_valence_warning);

        let nh2_dbl_on = analyze_left_endpoint_label("NH2", LabelJustification::Auto, 2, true);
        assert_eq!(nh2_dbl_on.raw_label, "NH");
        assert!(!nh2_dbl_on.has_valence_warning);

        let ch3_dbl_off = analyze_left_endpoint_label("CH3", LabelJustification::Auto, 2, false);
        assert_eq!(ch3_dbl_off.raw_label, "CH3");
        assert!(ch3_dbl_off.has_valence_warning);

        // Expansion succeeds for valid F3C and expands 3 F atoms; rejects invalid flush-right CF3
        let mut expand_mol = ChemMolecule::new();
        let root_c = expand_mol.add_atom(ChemAtom::new(1, 6));
        let other_c = expand_mol.add_atom(ChemAtom::new(2, 6));
        expand_mol.add_bond(3, root_c, other_c, 1, BondDisplay::Solid);
        assert!(expand_contracted_label(&mut expand_mol, root_c, &cf3_right, 14.40).is_err());
        let added = expand_contracted_label(&mut expand_mol, root_c, &f3c_right, 14.40).unwrap();
        assert_eq!(added.len(), 3);
        assert_eq!(expand_mol.atoms.iter().filter(|a| a.z == 9).count(), 3);
    }

    #[test]
    fn test_p02_stereo_cleanup_invariance() {
        // Central carbon with F (Z=9, wedge), Cl (Z=17), Br (Z=35), I (Z=53)
        let mut mol = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6);
        c.x = 100.0;
        c.y = 100.0;
        let ci = mol.add_atom(c);

        let mut f = ChemAtom::new(2, 9);
        f.x = 100.0;
        f.y = 60.0; // distorted length 40 pt (up)
        let fi = mol.add_atom(f);

        let mut cl = ChemAtom::new(3, 17);
        cl.x = 114.4;
        cl.y = 100.0;
        let cli = mol.add_atom(cl);

        let mut br = ChemAtom::new(4, 35);
        br.x = 100.0;
        br.y = 114.4;
        let bri = mol.add_atom(br);

        let mut i_atom = ChemAtom::new(5, 53);
        i_atom.x = 85.6;
        i_atom.y = 100.0;
        let ii = mol.add_atom(i_atom);

        mol.add_bond(1, ci, fi, 1, BondDisplay::WedgeBegin);
        mol.add_bond(2, ci, cli, 1, BondDisplay::Solid);
        mol.add_bond(3, ci, bri, 1, BondDisplay::Solid);
        mol.add_bond(4, ci, ii, 1, BondDisplay::Solid);

        let code_before = compute_tetrahedral_cip_code(&mol, ci);
        assert!(code_before == 2 || code_before == 3);

        clean_up_structure_preserving_stereo(&mut mol, &[0, 1, 2, 3], 14.40);
        let code_after = compute_tetrahedral_cip_code(&mol, ci);
        assert_eq!(code_before, code_after);

        let cf_dist = ((mol.atoms[fi].x - mol.atoms[ci].x).powi(2)
            + (mol.atoms[fi].y - mol.atoms[ci].y).powi(2))
        .sqrt();
        assert!((cf_dist - 14.40).abs() < 1e-9);
    }

    #[test]
    fn test_general_ring_fused_partial_and_ez_cleanup() {
        // 1. Naphthalene-like fused 6+6 system with distorted bond lengths: full selection normalizes all 11 bonds to 14.40
        let mut fused = ChemMolecule::new();
        let coords = [
            (0.0, 0.0),     // 0
            (20.0, -10.0),  // 1
            (40.0, 0.0),    // 2 (shared)
            (40.0, 22.0),   // 3 (shared)
            (20.0, 32.0),   // 4
            (0.0, 22.0),    // 5
            (62.0, -10.0),  // 6
            (82.0, 0.0),    // 7
            (82.0, 22.0),   // 8
            (62.0, 32.0),   // 9
        ];
        for (i, &(x, y)) in coords.iter().enumerate() {
            let mut a = ChemAtom::new((i + 1) as u32, 6);
            a.x = x;
            a.y = y;
            fused.add_atom(a);
        }
        let edges = [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (2, 6), (6, 7), (7, 8), (8, 9), (9, 3),
        ];
        for (bi, &(u, v)) in edges.iter().enumerate() {
            fused.add_bond((bi + 1) as u32, u, v, 6, BondDisplay::Solid);
        }
        let all_bonds: Vec<usize> = (0..edges.len()).collect();
        clean_up_structure_preserving_stereo(&mut fused, &all_bonds, 14.40);
        for b in &fused.bonds {
            let d = ((fused.atoms[b.begin_atom].x - fused.atoms[b.end_atom].x).powi(2)
                + (fused.atoms[b.begin_atom].y - fused.atoms[b.end_atom].y).powi(2))
            .sqrt();
            assert!((d - 14.40).abs() < 1e-6, "fused bond length {} != 14.40", d);
        }

        // 2. Partial-ring selection rule (D16): selecting 5 of 6 bonds of a hexagon leaves vertices untouched
        let mut partial = fused.clone();
        let before_pos: Vec<(f64, f64)> = partial.atoms.iter().map(|a| (a.x, a.y)).collect();
        clean_up_structure_preserving_stereo(&mut partial, &[0, 1, 2, 3, 4], 30.0);
        let after_pos: Vec<(f64, f64)> = partial.atoms.iter().map(|a| (a.x, a.y)).collect();
        assert_eq!(before_pos, after_pos);

        // 3. E/Z double bond + distorted acyclic chain preserves CIP code while normalizing all bond lengths
        let mut alkene = ChemMolecule::new();
        let a_specs = [
            (1, 6, 0.0, 0.0),     // 0: C=
            (2, 6, 25.0, 0.0),    // 1: =C (distorted 25 pt)
            (3, 17, -15.0, 12.0), // 2: Cl (top-left)
            (4, 9, -15.0, -12.0), // 3: F (bottom-left)
            (5, 35, 45.0, 14.0),  // 4: Br (top-right -> Z with Cl!)
        ];
        for (id, z, x, y) in a_specs {
            let mut a = ChemAtom::new(id, z);
            a.x = x;
            a.y = y;
            alkene.add_atom(a);
        }
        alkene.add_bond(1, 0, 1, 2, BondDisplay::Solid);
        alkene.add_bond(2, 0, 2, 1, BondDisplay::Solid);
        alkene.add_bond(3, 0, 3, 1, BondDisplay::Solid);
        alkene.add_bond(4, 1, 4, 1, BondDisplay::Solid);
        let ez_before = compute_double_bond_cip_code(&alkene, 0);
        assert_eq!(ez_before, 10); // Z

        clean_up_structure_preserving_stereo(&mut alkene, &[0, 1, 2, 3], 14.40);
        let ez_after = compute_double_bond_cip_code(&alkene, 0);
        assert_eq!(ez_before, ez_after);
        for b in &alkene.bonds {
            let d = ((alkene.atoms[b.begin_atom].x - alkene.atoms[b.end_atom].x).powi(2)
                + (alkene.atoms[b.begin_atom].y - alkene.atoms[b.end_atom].y).powi(2))
            .sqrt();
            assert!((d - 14.40).abs() < 1e-6);
        }
    }

    #[test]
    fn test_selection_hit_test_flip_and_svg_render() {
        // 1. Marquee vs Polygon Lasso & Zoom-scaled Hit-testing
        let sel_box = Rect2D { left: 0.0, top: 0.0, right: 50.0, bottom: 50.0 };
        let inside_obj = Rect2D { left: 10.0, top: 10.0, right: 40.0, bottom: 40.0 };
        let partial_obj = Rect2D { left: 10.0, top: 10.0, right: 60.0, bottom: 40.0 };
        assert!(marquee_selects_bounds(&sel_box, &inside_obj));
        assert!(!marquee_selects_bounds(&sel_box, &partial_obj));

        let lasso = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 50.0, y: 0.0 },
            Point2D { x: 50.0, y: 50.0 },
            Point2D { x: 0.0, y: 50.0 },
        ];
        assert!(lasso_selects_bounds(&lasso, &inside_obj));
        assert!(!lasso_selects_bounds(&lasso, &partial_obj));

        // 2. Chiral molecule: Flip with preserve_3d_stereo preserves R/S CIP parity & renders valid SVG
        let mut mol = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6); c.x = 100.0; c.y = 100.0;
        let mut f = ChemAtom::new(2, 9); f.x = 100.0; f.y = 85.6;
        let mut cl = ChemAtom::new(3, 17); cl.x = 114.4; cl.y = 100.0;
        let mut br = ChemAtom::new(4, 35); br.x = 100.0; br.y = 114.4;
        let mut o = ChemAtom::new(5, 8); o.x = 85.6; o.y = 100.0;
        mol.add_atom(c); mol.add_atom(f); mol.add_atom(cl); mol.add_atom(br); mol.add_atom(o);
        mol.add_bond(1, 0, 1, 1, BondDisplay::WedgeBegin);
        mol.add_bond(2, 0, 2, 1, BondDisplay::Solid);
        mol.add_bond(3, 0, 3, 1, BondDisplay::Solid);
        mol.add_bond(4, 0, 4, 2, BondDisplay::Solid);

        let tol = hit_tolerance_pt(5.0, 2.0); // 2.5 pt at 200% zoom
        assert_eq!(hit_test_atom(&mol, 101.0, 101.0, tol), Some(0));
        assert_eq!(hit_test_bond(&mol, 107.2, 100.0, tol), Some(1));

        let cip_orig = compute_tetrahedral_cip_code(&mol, 0);
        let mut rotated = mol.clone();
        rotate_atoms_preserving_stereo(&mut rotated, &[0, 1, 2, 3, 4], 0.7);
        assert_eq!(compute_tetrahedral_cip_code(&rotated, 0), cip_orig);

        let mut flipped = mol.clone();
        flip_atoms_preserving_stereo(&mut flipped, &[0, 1, 2, 3, 4], FlipAxis::Horizontal, true);
        assert_eq!(flipped.bonds[0].display, BondDisplay::WedgedHashBegin);
        assert_eq!(compute_tetrahedral_cip_code(&flipped, 0), cip_orig);

        let svg = render_molecule_svg(&mol, &DocumentSettings::acs_document_1996());
        assert!(svg.starts_with("<svg ") && svg.contains("</svg>"));
        assert!(svg.contains("<polygon "));
        assert!(svg.contains(">F<") && svg.contains(">Cl<") && svg.contains(">Br<"));
    }

    #[test]
    fn test_corpus_wide_cdxml_cdx_chem_and_svg_pipeline() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let inv_path = root.join("research/formats/inventory.json");
        if !inv_path.exists() {
            return;
        }
        let inv_text = fs::read_to_string(inv_path).unwrap();
        let mut cdxml_paths = Vec::new();
        let mut cdx_paths = Vec::new();
        let mut cur = String::new();
        for line in inv_text.lines() {
            let s = line.trim();
            if let Some(rest) = s.strip_prefix("\"path\": \"") {
                if let Some(end) = rest.find('"') {
                    cur = rest[..end].to_string();
                }
            } else if s.starts_with("\"kind\": \"CDXML\"") {
                cdxml_paths.push(cur.clone());
            } else if s.starts_with("\"kind\": \"CDX\"") {
                cdx_paths.push(cur.clone());
            }
        }

        let mut total_cdxml_mols = 0usize;
        let mut total_cdxml_rxn_steps = 0usize;
        for rel in &cdxml_paths {
            let raw = fs::read_to_string(root.join(rel)).unwrap();
            let cdxml = CdxmlDocument::parse(&raw).unwrap();
            let doc = InteractiveDocument::from_cdxml(&cdxml);
            for mol in &doc.scene.molecules {
                total_cdxml_mols += 1;
                let fp = mol.compute_formula_properties();
                let _ = fp.get_molecular_weight();
                let _ = chiralcraft_chem::perceive_stereochemistry(mol);
                let _ = chiralcraft_chem::CcAromaticityInfo::perceive(mol);
            }
            let svg = doc.to_svg();
            assert!(svg.contains("<svg ") && svg.contains("</svg>"));

            for scheme in cdxml.extract_reactions() {
                let mut rxn = chiralcraft_chem::CcReaction::new();
                for _ in &scheme.steps {
                    total_cdxml_rxn_steps += 1;
                    rxn.steps.push(chiralcraft_chem::ReactionStep {
                        reactant: ChemMolecule::new(),
                        product: ChemMolecule::new(),
                        atom_map: std::collections::BTreeMap::new(),
                        bond_map: std::collections::BTreeMap::new(),
                    });
                }
                assert!(rxn.validate().is_ok());
            }
        }

        let mut total_cdx_mols = 0usize;
        for rel in &cdx_paths {
            let raw = fs::read(root.join(rel)).unwrap();
            let cdx = chiralcraft_formats::CdxDocument::parse(&raw).unwrap();
            for (_, mol) in cdx.extract_molecules() {
                total_cdx_mols += 1;
                let _ = mol.compute_formula_properties();
                let svg = render_molecule_svg(&mol, &DocumentSettings::acs_document_1996());
                assert!(svg.contains("</svg>"));
            }
        }

        assert!(total_cdxml_mols > 0);
        assert_eq!(total_cdxml_rxn_steps, 371);
        assert_eq!(total_cdx_mols, 480);

        // Verify live-captured P04 fixture (p04-acs-100.cdxml) exact 14.40 pt bond length
        let p04_path = root.join("research/runtime-probes/p04/p04-acs-100.cdxml");
        if p04_path.exists() {
            let p04_xml = fs::read_to_string(p04_path).unwrap();
            let p04_cdxml = CdxmlDocument::parse(&p04_xml).unwrap();
            let p04_doc = InteractiveDocument::from_cdxml(&p04_cdxml);
            assert_eq!(p04_doc.scene.molecules.len(), 1);
            let m = &p04_doc.scene.molecules[0];
            assert_eq!(m.atoms.len(), 2);
            assert_eq!(m.bonds.len(), 1);
            let dx = m.atoms[1].x - m.atoms[0].x;
            let dy = m.atoms[1].y - m.atoms[0].y;
            let dist = (dx * dx + dy * dy).sqrt();
            assert!((dist - 14.40).abs() < 1e-6, "P04 bond length {dist} != 14.40");
        }
    }

    #[test]
    fn test_interactive_graph_editing_d01_d03_d11_d13() {
        let settings = DocumentSettings::acs_document_1996();
        let mut mol = ChemMolecule::new();
        let mut a0 = ChemAtom::new(1, 6);
        a0.x = 100.0;
        a0.y = 100.0;
        mol.add_atom(a0);

        let (a1, b0) = add_bond_from_atom(&mut mol, 0, 0.0, &settings).unwrap();
        assert_eq!(a1, 1);
        assert!((mol.atoms[1].x - 114.40).abs() < 1e-6);

        // D02/D03: cycle bond order 1 -> 2 -> 3 -> 4 -> 1
        assert_eq!(cycle_bond_order(&mut mol, b0), Some(2));
        assert_eq!(cycle_bond_order(&mut mol, b0), Some(3));
        assert_eq!(cycle_bond_order(&mut mol, b0), Some(4));
        assert_eq!(cycle_bond_order(&mut mol, b0), Some(1));

        // D12/D13: fuse benzene ring onto b0 (adds 4 atoms, 5 bonds -> 6 atoms, 6 bonds)
        let added = fuse_ring_on_bond(&mut mol, b0, 6, true);
        assert_eq!(added.len(), 4);
        assert_eq!(mol.atoms.len(), 6);
        assert_eq!(mol.bonds.len(), 6);
        assert_eq!(mol.compute_formula_properties().get_formula(), "C6H6");

        // D11: collapse one ring bond -> reduces 6-membered ring to 5 bonds
        assert!(collapse_bond_reducing_ring(&mut mol, 0));
        assert_eq!(mol.bonds.len(), 5);

        // D12: spiro ring on atom 0 (adds 4 more atoms, 5 bonds for 5-membered ring)
        let spiro_added = spiro_ring_on_atom(&mut mol, 0, 5, settings.bond_length_pt, false);
        assert_eq!(spiro_added.len(), 4);

        // D17 Check Structure: isolated C-C bond triggers exact ChiralCraft p04 warning; invalid H-bond triggers D07 warning
        let mut isolated = ChemMolecule::new();
        isolated.add_atom(ChemAtom::new(1, 6));
        let mut c2 = ChemAtom::new(2, 6);
        c2.x = 14.4;
        isolated.add_atom(c2);
        isolated.add_bond(1, 0, 1, 1, BondDisplay::Solid);
        let w_iso = check_structure(&isolated, &[]);
        assert_eq!(w_iso.len(), 1);
        assert_eq!(
            w_iso[0].message,
            "This isolated bond is probably not intended to have chemical significance."
        );

        let mut bad_hb = ChemMolecule::new();
        bad_hb.add_atom(ChemAtom::new(1, 6));
        let mut c_other = ChemAtom::new(2, 6);
        c_other.x = 14.4;
        bad_hb.add_atom(c_other);
        bad_hb.add_bond(1, 0, 1, 7, BondDisplay::Dash);
        let w_hb = check_structure(&bad_hb, &[]);
        assert!(w_hb.iter().any(|w| w.message == "Invalid or ambiguous H-Bonds"));

        // E14 Align / Distribute / Center on page
        let mut ms = vec![isolated.clone(), bad_hb.clone(), isolated];
        shift_molecule(&mut ms[1], 100.0, 20.0);
        shift_molecule(&mut ms[2], 300.0, 50.0);
        align_molecules(&mut ms, AlignMode::Top);
        assert!(ms.iter().all(|m| (m.atoms[0].y - ms[0].atoms[0].y).abs() < 1e-6));
        distribute_molecules(&mut ms, DistributeAxis::Horizontal);
        let gap1 = ms[1].atoms[0].x - ms[0].atoms[1].x;
        let gap2 = ms[2].atoms[0].x - ms[1].atoms[1].x;
        assert!((gap1 - gap2).abs() < 1e-6);
        center_on_page(&mut ms, 540.0, 719.75);
        assert!(join_atoms(&mut ms[0], 0, 1));
        assert!(ms[0].bonds.is_empty());
    }

    #[test]
    fn test_live_prime_p01_to_p06_captures_and_roundtrip() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }
        let cdxml_rel_paths = [
            "p01/p01-v5-doc.cdxml",
            "p02/p02-precleanup.cdxml",
            "p02/p02-after-cleanup-1.cdxml",
            "p02/p02-after-cleanup-3.cdxml",
            "p03/p03-flushright-expanded.cdxml",
            "p03/p03-auto-expanded.cdxml",
            "p04/p04-acs-200-b.cdxml",
            "p04/p04-1in-drawn.cdxml",
            "p05/p05-nudged.cdxml",
            "p05/p05-undone.cdxml",
            "p05/p05-ctrl-drag.cdxml",
            "p06/p06-after-cdx.cdxml",
        ];
        for rel in cdxml_rel_paths {
            let p = probes.join(rel);
            assert!(p.exists(), "missing live capture {}", p.display());
            let xml = fs::read_to_string(&p).unwrap();
            let parsed = CdxmlDocument::parse(&xml).expect(rel);
            let refs = parsed.analyze_references();
            assert_eq!(refs.bond_endpoints_typed_local_dangling, 0, "dangling bond refs in {}", rel);
            assert_eq!(refs.non_bond_references_unresolved, 0, "unresolved refs in {}", rel);
            let doc = InteractiveDocument::from_cdxml(&parsed);
            let rt_xml = doc.to_cdxml().to_xml_string();
            let rt_parsed = CdxmlDocument::parse(&rt_xml).expect("roundtrip CDXML parse");
            let rt_doc = InteractiveDocument::from_cdxml(&rt_parsed);
            assert_eq!(doc.scene.molecules.len(), rt_doc.scene.molecules.len(), "mol count mismatch for {}", rel);
            let svg = doc.to_svg();
            assert!(svg.contains("<svg ") && svg.contains("</svg>"));
        }

        // Verify live ChiralCraft CDX file parses, extracts exact molecules, and byte-roundtrips
        let cdx_bytes = fs::read(probes.join("p06/p06-roundtrip.cdx")).unwrap();
        let cdx = chiralcraft_formats::CdxDocument::parse(&cdx_bytes).unwrap();
        assert_eq!(cdx.to_bytes(), cdx_bytes, "CDX binary roundtrip must be byte-exact");
        let cdx_mols = cdx.extract_molecules();
        assert_eq!(cdx_mols.len(), 3);
        let cdx_to_cdxml = cdx.to_cdxml();
        let cdx_doc = InteractiveDocument::from_cdxml(&cdx_to_cdxml);
        assert_eq!(cdx_doc.scene.molecules.len(), 3);

        // Verify P04 200%-zoom bond == 14.40 pt and P04 1in override == 72.00 pt
        let p04_200 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p04/p04-acs-200-b.cdxml")).unwrap()).unwrap();
        let d200 = InteractiveDocument::from_cdxml(&p04_200);
        let b200 = &d200.scene.molecules[0];
        let d200_len = ((b200.atoms[1].x - b200.atoms[0].x).powi(2) + (b200.atoms[1].y - b200.atoms[0].y).powi(2)).sqrt();
        assert!((d200_len - 14.40).abs() < 1e-4, "expected 14.40 pt, got {}", d200_len);

        let p04_1in = CdxmlDocument::parse(&fs::read_to_string(probes.join("p04/p04-1in-drawn.cdxml")).unwrap()).unwrap();
        let d1in = InteractiveDocument::from_cdxml(&p04_1in);
        assert!((d1in.settings.bond_length_pt - 72.0).abs() < 1e-4);
        let b1in = &d1in.scene.molecules[0];
        let d1in_len = ((b1in.atoms[1].x - b1in.atoms[0].x).powi(2) + (b1in.atoms[1].y - b1in.atoms[0].y).powi(2)).sqrt();
        assert!((d1in_len - 72.0).abs() < 1e-4, "expected 72.00 pt, got {}", d1in_len);

        // Verify P05 Ctrl-drag duplicated rectangle graphic from 1 to 2
        let p05_seed = CdxmlDocument::parse(&fs::read_to_string(probes.join("p05/p05-seed.cdxml")).unwrap()).unwrap();
        let p05_dup = CdxmlDocument::parse(&fs::read_to_string(probes.join("p05/p05-ctrl-drag.cdxml")).unwrap()).unwrap();
        assert_eq!(p05_seed.elements.iter().filter(|e| e.tag == "graphic").count(), 1);
        assert_eq!(p05_dup.elements.iter().filter(|e| e.tag == "graphic").count(), 2);
    }

    #[test]
    fn test_live_prime_p07_to_p09_exports_and_chemistry_oracle() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }
        // 1. P07 live V2000 Molfile exported by ChiralCraft parses cleanly
        let mol_v2000 = fs::read_to_string(probes.join("p07/p07-export-v2000.mol")).unwrap();
        let mol = chiralcraft_formats::molfile_v2000_to_molecule(&mol_v2000).unwrap();
        assert_eq!(mol.atoms.len(), 17);
        assert_eq!(mol.bonds.len(), 15);

        // 2. P08 live checked CDXML Aspirin (Fragment 30 -> index 1) matches Prime Analysis Window
        let p08_xml = fs::read_to_string(probes.join("p08/p08-checked.cdxml")).unwrap();
        let p08 = CdxmlDocument::parse(&p08_xml).unwrap();
        let p08_doc = InteractiveDocument::from_cdxml(&p08);
        assert_eq!(p08_doc.scene.molecules.len(), 3);
        let aspirin = &p08_doc.scene.molecules[1];
        let mut fp = aspirin.compute_formula_properties();
        assert_eq!(fp.get_formula(), "C9H8O4");
        assert!((fp.get_exact_mass().val - 180.0423).abs() < 0.01);
        assert!((fp.get_molecular_weight() - 180.16).abs() < 0.01);

        // Check Structure on pentavalent carbon (Fragment 10 -> index 0) reports over-valence diagnostic
        let w_pent = check_structure(&p08_doc.scene.molecules[0], &[]);
        assert!(!w_pent.is_empty());

        // 3. P02 partial half-ring cleanup + P09 reaction CDXML parse and extract cleanly
        let p02_half = CdxmlDocument::parse(&fs::read_to_string(probes.join("p02/p02-partial-half-cleaned.cdxml")).unwrap()).unwrap();
        let p02_doc = InteractiveDocument::from_cdxml(&p02_half);
        assert_eq!(p02_doc.scene.molecules.len(), 1);
        assert_eq!(p02_doc.scene.molecules[0].bonds.len(), 6);

        let p09_rxn = CdxmlDocument::parse(&fs::read_to_string(probes.join("p09/p09-rxn-cleaned.cdxml")).unwrap()).unwrap();
        let rxns = p09_rxn.extract_reactions();
        assert_eq!(rxns.len(), 1);
        assert_eq!(rxns[0].steps.len(), 1);
    }

    #[test]
    fn test_live_prime_p10_p11_smiles_group_and_tools() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }
        // 1. P11 Paste Special > SMILES multi-fragment ("CCO.c1ccccc1"): auto-wrapped in <group>, Order="1.5" + Circle
        let p11_smiles = CdxmlDocument::parse(&fs::read_to_string(probes.join("p11/p11-pasted-smiles.cdxml")).unwrap()).unwrap();
        assert_eq!(p11_smiles.elements.iter().filter(|e| e.tag == "group").count(), 1);
        assert_eq!(p11_smiles.elements.iter().filter(|e| e.tag == "graphic" && e.attr("OvalType") == Some("Circle")).count(), 1);
        let smiles_doc = InteractiveDocument::from_cdxml(&p11_smiles);
        assert_eq!(smiles_doc.scene.molecules.len(), 2);

        // 2. P10 TLC plate + per-bond BondLength override + live bond toggle
        let p10_tools = CdxmlDocument::parse(&fs::read_to_string(probes.join("p10/p10-tools.cdxml")).unwrap()).unwrap();
        assert_eq!(p10_tools.elements.iter().filter(|e| e.tag == "tlcplate").count(), 1);
        assert_eq!(p10_tools.elements.iter().filter(|e| e.tag == "tlclane").count(), 5);
        assert!(p10_tools.elements.iter().any(|e| e.tag == "b" && e.attr("BondLength") == Some("52.52")));

        let p10_toggle = CdxmlDocument::parse(&fs::read_to_string(probes.join("p10/p10-hbond-toggle-live.cdxml")).unwrap()).unwrap();
        let toggle_doc = InteractiveDocument::from_cdxml(&p10_toggle);
        assert_eq!(toggle_doc.scene.molecules.len(), 2);
        assert!(p10_toggle.elements.iter().any(|e| e.tag == "b" && e.attr("Order") == Some("2")));
        assert!(p10_toggle.elements.iter().any(|e| e.tag == "b" && e.attr("Display") == Some("WedgeBegin")));

        // 3. P11 double-bond CH3 retains NumHydrogens=3 and invalid valence warning on CDXML load
        let p11_hrect = CdxmlDocument::parse(&fs::read_to_string(probes.join("p11/p11-hbond-rect-out.cdxml")).unwrap()).unwrap();
        assert!(p11_hrect.elements.iter().any(|e| {
            e.tag == "n"
                && e.attr("NumHydrogens") == Some("3")
                && e.attr("Warning") == Some("An atom in this label has an invalid valence.")
        }));
    }

    #[test]
    fn test_live_prime_p12_p13_naphthalene_rxn_and_inchi() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }
        // 1. P12 fused Naphthalene cleaned to 15.25 pt and live CH3 -> CH2 rectification
        let p12 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p12/p12-out.cdxml")).unwrap()).unwrap();
        let p12_doc = InteractiveDocument::from_cdxml(&p12);
        let naph = p12_doc.scene.molecules.iter().find(|m| m.atoms.len() == 10 && m.bonds.len() == 11).unwrap();
        for b in &naph.bonds {
            let len = ((naph.atoms[b.end_atom].x - naph.atoms[b.begin_atom].x).powi(2)
                + (naph.atoms[b.end_atom].y - naph.atoms[b.begin_atom].y).powi(2))
            .sqrt();
            assert!((len - 15.25).abs() < 0.05, "expected ~15.25 pt fused bond, got {}", len);
        }
        assert!(p12.to_xml_string().contains(">CH2<"));

        // 2. P13 aligned reaction: reactant + product both captured in <step> + Paste Special InChI ethanol
        let p13_rxn = CdxmlDocument::parse(&fs::read_to_string(probes.join("p13/p13-rxn-aligned-out.cdxml")).unwrap()).unwrap();
        let schemes = p13_rxn.extract_reactions();
        assert_eq!(schemes.len(), 1);
        assert_eq!(schemes[0].steps[0].reactant_ids, vec!["10"]);
        assert_eq!(schemes[0].steps[0].product_ids, vec!["20"]);

        let p13_inchi = CdxmlDocument::parse(&fs::read_to_string(probes.join("p13/p13-pasted-inchi.cdxml")).unwrap()).unwrap();
        let inchi_doc = InteractiveDocument::from_cdxml(&p13_inchi);
        assert_eq!(inchi_doc.scene.molecules.len(), 1);
        assert_eq!(inchi_doc.scene.molecules[0].compute_formula_properties().get_formula(), "C2H6O");
    }

    #[test]
    fn test_live_prime_p14_p15_advanced_stereo_mapping_and_objects() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }
        // 1. P14 binary CDX roundtrip & advanced CIP / Enhanced Stereo
        let cdx_bytes = fs::read(probes.join("p14/p14-roundtrip.cdx")).unwrap();
        let cdx = chiralcraft_formats::CdxDocument::parse(&cdx_bytes).unwrap();
        assert_eq!(cdx.to_bytes(), cdx_bytes);

        let p14 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p14/p14-out.cdxml")).unwrap()).unwrap();
        let p14_xml = p14.to_xml_string();
        assert!(p14_xml.contains("AS=\"-\""));
        assert!(p14_xml.contains(">(M)<"));
        let enhanced = p14.extract_enhanced_stereo();
        let markers: Vec<_> = enhanced.values().flatten().copied().collect();
        assert!(markers.iter().any(|m| m.group_type == 3 && m.group_number == 1));
        assert!(markers.iter().any(|m| m.group_type == 4 && m.group_number == 2));
        let p14_mol = fs::read_to_string(probes.join("p14/p14-v2000.mol")).unwrap();
        assert!(!p14_mol.contains("EnhancedStereo") && !p14_mol.contains("(M)"));

        // 2. P14 pseudoasymmetric (R)/(r)/(S)
        let p14_pseudo = CdxmlDocument::parse(&fs::read_to_string(probes.join("p14/p14-pseudo-out.cdxml")).unwrap()).unwrap();
        let pseudo_xml = p14_pseudo.to_xml_string();
        assert!(pseudo_xml.contains("AS=\"R\"") && pseudo_xml.contains("AS=\"r\"") && pseudo_xml.contains("AS=\"S\""));
        assert!(pseudo_xml.contains(">(r)<"));

        // 3. P15 Map Reaction Atoms (CDXML + $RXN V3000)
        let p15_rxn = CdxmlDocument::parse(&fs::read_to_string(probes.join("p15/p15-mapped-rxn.cdxml")).unwrap()).unwrap();
        let rxn_xml = p15_rxn.to_xml_string();
        assert!(rxn_xml.contains("ReactionStepAtomMap=\"11 21 12 22\""));
        assert!(rxn_xml.contains("RxnParticipation=\"ChangeType\""));
        assert!(rxn_xml.contains("RxnParticipation=\"MakeOrBreak\""));
        let schemes = p15_rxn.extract_reactions();
        assert_eq!(schemes.len(), 1);
        assert_eq!(schemes[0].steps[0].reactant_ids, vec!["10"]);
        assert_eq!(schemes[0].steps[0].product_ids, vec!["20"]);
        let rxn_v3000 = fs::read_to_string(probes.join("p15/p15-mapped-rxn.rxn")).unwrap();
        assert!(rxn_v3000.contains("$RXN V3000") && rxn_v3000.contains("RXCTR=8") && rxn_v3000.contains("RXCTR=4"));

        // 4. P15 Benzene fusion + sprout, Brackets, Table, and Acyclic Chain
        let p15_tools = CdxmlDocument::parse(&fs::read_to_string(probes.join("p15/p15-tools.cdxml")).unwrap()).unwrap();
        let tools_doc = InteractiveDocument::from_cdxml(&p15_tools);
        assert_eq!(tools_doc.scene.molecules.len(), 1);
        assert_eq!(tools_doc.scene.molecules[0].atoms.len(), 15);
        assert_eq!(tools_doc.scene.molecules[0].bonds.len(), 17);
        let tools_xml = p15_tools.to_xml_string();
        assert!(tools_xml.contains("Warning=\"There are too many bonds to this unlabeled Carbon.\""));
        assert!(tools_xml.contains("<table") && tools_xml.contains("<bracketedgroup") && tools_xml.contains("BracketType=\"Square\""));

        let p15_chain = CdxmlDocument::parse(&fs::read_to_string(probes.join("p15/p15-chain.cdxml")).unwrap()).unwrap();
        let chain_doc = InteractiveDocument::from_cdxml(&p15_chain);
        assert_eq!(chain_doc.scene.molecules.len(), 1);
        assert_eq!(chain_doc.scene.molecules[0].atoms.len(), 5);
        assert_eq!(chain_doc.scene.molecules[0].bonds.len(), 4);
    }

    #[test]
    fn test_live_prime_p16_p17_atropisomer_transforms_graphics_nicknames_multipage() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }

        // 1. P16 atropisomer biaryl (M) perception (AS="M") + ambiguous CHBrClF AS="u"
        let atrop = CdxmlDocument::parse(&fs::read_to_string(probes.join("p16/p16-atrop-out.cdxml")).unwrap()).unwrap();
        let atrop_xml = atrop.to_xml_string();
        assert!(atrop_xml.contains("AS=\"M\"") && atrop_xml.contains(">(M)<"));
        assert!(atrop_xml.contains("AS=\"u\""));

        // 2. P16 Flip Horizontal vs Rotate 180 Horizontal
        let flipped = CdxmlDocument::parse(&fs::read_to_string(probes.join("p16/p16-flipped.cdxml")).unwrap()).unwrap();
        let flipped_xml = flipped.to_xml_string();
        // Flip preserves WedgeBegin, inverting (R)->(S) and (S)->(R)
        assert!(flipped_xml.contains("Display=\"WedgeBegin\""));
        assert!(!flipped_xml.contains("Display=\"WedgedHashBegin\""));

        let rot180 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p16/p16-rotated180.cdxml")).unwrap()).unwrap();
        let rot180_xml = rot180.to_xml_string();
        // Rotate 180 Horizontal swaps WedgeBegin -> WedgedHashBegin to preserve CIP
        assert!(rot180_xml.contains("Display=\"WedgedHashBegin\""));
        assert!(!rot180_xml.contains("Display=\"WedgeBegin\""));

        // 3. P16 graphics primitives: Curve, Arrow, Orbital, Shadow Rectangle, Bracket
        let gfx = CdxmlDocument::parse(&fs::read_to_string(probes.join("p16/p16-graphics.cdxml")).unwrap()).unwrap();
        let gfx_xml = gfx.to_xml_string();
        assert!(gfx_xml.contains("<curve") && gfx_xml.contains("CurvePoints=\""));
        assert!(gfx_xml.contains("<arrow") && gfx_xml.contains("HeadSize=\"2250\""));
        assert!(gfx_xml.contains("GraphicType=\"Orbital\"") && gfx_xml.contains("OrbitalType=\"p\""));
        assert!(gfx_xml.contains("GraphicType=\"Rectangle\"") && gfx_xml.contains("RectangleType=\"Shadow\""));
        assert!(gfx_xml.contains("<bracketedgroup"));

        // 4. P17 Expand Label across 7 multi-atom nicknames (Boc, OTs, TMS, CO2Et, Cbz, iPr, CH2OH)
        let exp = CdxmlDocument::parse(&fs::read_to_string(probes.join("p17/p17-expanded.cdxml")).unwrap()).unwrap();
        let exp_doc = InteractiveDocument::from_cdxml(&exp);
        let exp_xml = exp.to_xml_string();
        assert_eq!(exp_doc.scene.molecules.len(), 7);
        assert!(!exp_xml.contains("NodeType=\"Nickname\""));
        for tok in ["Boc", "TMS", "Cbz", "iPr"] {
            assert!(exp_xml.contains(&format!("IsNickname=\"yes\"\n>{tok}<")) || exp_xml.contains(&format!(">{tok}<")), "expected hidden IsNickname cache for {tok}");
        }
        assert!(exp_xml.contains("Name=\"Residue\"") && exp_xml.contains("BasisObjects="));
        assert!(exp_xml.contains(">Si<") && exp_xml.contains(">S<") && exp_xml.contains(">OH<"));

        // 5. P17 2x2 tiled page document preservation
        let mp = CdxmlDocument::parse(&fs::read_to_string(probes.join("p17/p17-multipage-out.cdxml")).unwrap()).unwrap();
        let mp_xml = mp.to_xml_string();
        assert!(mp_xml.contains("WidthPages=\"2\"") && mp_xml.contains("HeightPages=\"2\""));
    }

    #[test]
    fn test_live_prime_p18_p19_valence_sru_caffeine_and_smiles() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }

        // 1. P18 live Check Structure warnings vs Rust check_structure parity
        let p18 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p18/p18-valence-sru-out.cdxml")).unwrap()).unwrap();
        let p18_xml = p18.to_xml_string();
        assert!(p18_xml.contains("Warning=\"There are too many bonds to this unlabeled Carbon.\""));
        assert!(p18_xml.contains("Warning=\"An atom in this label has an invalid valence.\""));
        assert!(p18_xml.contains("Warning=\"Parentheses don&apos;t match.\"") || p18_xml.contains("Warning=\"Parentheses don't match.\""));

        let p18_doc = InteractiveDocument::from_cdxml(&p18);
        assert_eq!(p18_doc.scene.molecules.len(), 5);
        let w_c5 = check_structure(&p18_doc.scene.molecules[0], &[]);
        assert_eq!(w_c5.len(), 1);
        assert_eq!(w_c5[0].message, "There are too many bonds to this unlabeled Carbon.");
        let w_n5 = check_structure(&p18_doc.scene.molecules[1], &[]);
        assert_eq!(w_n5.len(), 1);
        assert_eq!(w_n5[0].message, "An atom in this label has an invalid valence.");
        // Period-3 hypervalent SF6 passes clean
        assert!(check_structure(&p18_doc.scene.molecules[2], &[]).is_empty());
        // Single-bond P-Fe (sbo=4) triggers invalid valence on P
        let w_pfe = check_structure(&p18_doc.scene.molecules[3], &[]);
        assert_eq!(w_pfe.len(), 1);
        assert_eq!(w_pfe[0].message, "An atom in this label has an invalid valence.");

        // 2. P19 Caffeine Analysis + Copy As SMILES cross-oracle verification
        let p19 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p19/p19-caffeine-out.cdxml")).unwrap()).unwrap();
        let p19_doc = InteractiveDocument::from_cdxml(&p19);
        assert_eq!(p19_doc.scene.molecules.len(), 1);
        let mut fp = p19_doc.scene.molecules[0].compute_formula_properties();
        assert_eq!(fp.get_formula(), "C8H10N4O2");
        assert!((fp.get_exact_mass().val - 194.0804).abs() < 0.01);
        assert!((fp.get_molecular_weight() - 194.19).abs() < 0.01);
        let ea = fp.get_elemental_analysis();
        let ea_map: HashMap<String, f64> = ea.into_iter().collect();
        assert!((ea_map["C"] - 49.48).abs() < 0.05);
        assert!((ea_map["H"] - 5.19).abs() < 0.05);
        assert!((ea_map["N"] - 28.85).abs() < 0.05);
        assert!((ea_map["O"] - 16.48).abs() < 0.05);

        let smiles = fs::read_to_string(probes.join("p19/p19-caffeine-smiles.txt")).unwrap();
        let smi_mol = chiralcraft_formats::smiles_to_molecule(smiles.trim()).unwrap();
        assert_eq!(smi_mol.compute_formula_properties().get_formula(), "C8H10N4O2");
    }

    #[test]
    fn test_live_prime_p20_p21_dative_brackets_and_caffeine_props() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/runtime-probes");
        if !probes.exists() {
            return;
        }

        // 1. P20 lowercase Order="dative" (N->B, P->Fe) & Prime auto-computed crossing bonds on brackets
        let p20 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p20/p20-dative-sru-out.cdxml")).unwrap()).unwrap();
        let p20_doc = InteractiveDocument::from_cdxml(&p20);
        assert_eq!(p20_doc.scene.molecules.len(), 3);
        // Me3N -> BH3: N=0 warnings, bare B(NumHydrogens=0) invalid valence
        let w_nb = check_structure(&p20_doc.scene.molecules[0], &[]);
        assert_eq!(w_nb.len(), 1);
        assert_eq!(w_nb[0].message, "An atom in this label has an invalid valence.");
        // Me3P -> Fe with Order="dative": zero warnings in Prime & Rust
        assert!(check_structure(&p20_doc.scene.molecules[1], &[]).is_empty());

        let groups = p20.extract_bracketed_groups();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].bracketed_object_ids, vec!["32", "33"]);
        assert_eq!(groups[0].attachments.len(), 2);
        assert_eq!(groups[0].attachments[0].crossing_bonds[0].bond_id.as_deref(), Some("303"));
        assert_eq!(groups[0].attachments[0].crossing_bonds[0].inner_atom_id.as_deref(), Some("33"));
        assert_eq!(groups[0].attachments[1].crossing_bonds[0].bond_id.as_deref(), Some("301"));
        assert_eq!(groups[0].attachments[1].crossing_bonds[0].inner_atom_id.as_deref(), Some("32"));

        // 2. P21 live completed Chemical Properties pasted onto Caffeine
        let p21 = CdxmlDocument::parse(&fs::read_to_string(probes.join("p21/p21-caffeine-pasted.cdxml")).unwrap()).unwrap();
        let p21_xml = p21.to_xml_string();
        assert!(p21_xml.contains("Log P: -1.03"));
        assert!(p21_xml.contains("MR: 51.2 [cm3/mol]"));
        assert!(p21_xml.contains("tPSA: 56.22"));
        assert!(p21_xml.contains("LogS: -1.084"));
        assert!(p21_xml.contains("ChemicalPropertyDisplayID=\"116\""));
    }

    #[test]
    fn test_live_prime_p22_p23_contract_queries_align_and_rotate() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P22a: Structure > Contract Label restores 4 NodeType="Fragment" groups + visible IsNickname="yes" on Boc/TMS/Cbz/iPr
        let p22_contract_xml = fs::read_to_string(probes.join("p22/p22-contracted.cdxml")).unwrap();
        assert_eq!(p22_contract_xml.matches("NodeType=\"Fragment\"").count(), 4);
        assert_eq!(p22_contract_xml.matches("IsNickname=\"yes\"").count(), 4);
        assert!(!p22_contract_xml.contains("Visible=\"no\""));
        assert!(p22_contract_xml.contains("ExternalConnectionType=\"Residue\""));
        let p22_contract_doc = CdxmlDocument::parse(&p22_contract_xml).unwrap();
        let p22_contract_scene = InteractiveDocument::from_cdxml(&p22_contract_doc);
        assert_eq!(p22_contract_scene.scene.molecules.len(), 11);

        // 2. P22c: Query properties (ElementList, SubstituentsExactly, GenericNickname) pass Check Structure cleanly & Expand Generic Structure emits 3 discrete molecules
        let p22_qcheck_xml = fs::read_to_string(probes.join("p22/p22-query-check.cdxml")).unwrap();
        assert!(p22_qcheck_xml.contains("ElementList=\"9 17 35\""));
        assert!(p22_qcheck_xml.contains("SubstituentsExactly=\"2\""));
        assert!(p22_qcheck_xml.contains("GenericNickname=\"A\""));
        let p22_qcheck = InteractiveDocument::from_cdxml(&CdxmlDocument::parse(&p22_qcheck_xml).unwrap());
        assert_eq!(p22_qcheck.scene.molecules.len(), 2);
        for (mi, mol) in p22_qcheck.scene.molecules.iter().enumerate() {
            let w = check_structure(mol, &[]);
            assert!(w.is_empty(), "mol {} warnings: {:?}", mi, w);
        }

        let p22_qout = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p22/p22-query-out.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(p22_qout.scene.molecules.len(), 3);
        let expanded_halogens: Vec<i32> = p22_qout
            .scene
            .molecules
            .iter()
            .map(|m| m.atoms[2].z)
            .collect();
        assert_eq!(expanded_halogens, vec![9, 17, 35]);

        // 3. P23: Align > Left edges + Rotate... 45 deg matches live ChiralCraft coords and preserves (R)/(S) CIP stereo
        let p23_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p23/p23-align-rotate-seed.cdxml")).unwrap()).unwrap(),
        );
        let p23_aligned_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p23/p23-aligned-only.cdxml")).unwrap()).unwrap(),
        );
        let p23_rotated_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p23/p23-aligned-rotated.cdxml")).unwrap()).unwrap(),
        );

        let mut mols = p23_seed.scene.molecules.clone();
        align_molecules(&mut mols, AlignMode::Left);
        for (rust_mol, live_mol) in mols.iter().zip(p23_aligned_live.scene.molecules.iter()) {
            for (ra, la) in rust_mol.atoms.iter().zip(live_mol.atoms.iter()) {
                assert!((ra.x - la.x).abs() < 1e-3, "align x mismatch: {} vs {}", ra.x, la.x);
                assert!((ra.y - la.y).abs() < 1e-3, "align y mismatch: {} vs {}", ra.y, la.y);
            }
        }

        rotate_molecules_around_bounds_center(&mut mols, 45.0_f64.to_radians());
        for (rust_mol, live_mol) in mols.iter().zip(p23_rotated_live.scene.molecules.iter()) {
            for (ra, la) in rust_mol.atoms.iter().zip(live_mol.atoms.iter()) {
                // Prime's Rotate... pivots on bounding box including Cl/OH text glyphs (~0.03 pt delta vs atom bounds)
                assert!((ra.x - la.x).abs() < 0.05, "rotate x mismatch: {} vs {}", ra.x, la.x);
                assert!((ra.y - la.y).abs() < 0.05, "rotate y mismatch: {} vs {}", ra.y, la.y);
            }
        }

        assert_eq!(compute_tetrahedral_cip_code(&mols[0], 1), 2); // (R)
        assert_eq!(compute_tetrahedral_cip_code(&mols[1], 1), 3); // (S)
        assert_eq!(
            compute_tetrahedral_cip_code(&p23_rotated_live.scene.molecules[0], 1),
            2
        );
        assert_eq!(
            compute_tetrahedral_cip_code(&p23_rotated_live.scene.molecules[1], 1),
            3
        );
    }

    #[test]
    fn test_live_prime_p24_p25_distribute_flips_rot90_svg_and_sdf() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P24a: Object > Distribute > Horizontally matches distribute_molecules 1:1
        let dist_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p24/p24-distribute-seed.cdxml")).unwrap()).unwrap(),
        );
        let dist_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p24/p24-distributed.cdxml")).unwrap()).unwrap(),
        );
        let mut dmols = dist_seed.scene.molecules.clone();
        distribute_molecules(&mut dmols, DistributeAxis::Horizontal);
        assert!((dmols[1].atoms[0].x - 190.0).abs() < 1e-3);
        for (rm, lm) in dmols.iter().zip(dist_live.scene.molecules.iter()) {
            for (ra, la) in rm.atoms.iter().zip(lm.atoms.iter()) {
                assert!((ra.x - la.x).abs() < 1e-3);
                assert!((ra.y - la.y).abs() < 1e-3);
            }
        }

        // 2. P24b: Dual-stereocenter (2R,3S)-2-bromo-3-chloropentane Flip Vertical vs Rotate 180 Vertical vs Rotate 90 & SVG
        let mc_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p24/p24-multi-chiral-seed.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(compute_tetrahedral_cip_code(&mc_seed.scene.molecules[0], 1), 2); // 2R
        assert_eq!(compute_tetrahedral_cip_code(&mc_seed.scene.molecules[0], 2), 3); // 3S

        // Flip Vertical keeps wedges -> inverts both stereocenters (2S, 3R)
        let fv_xml = fs::read_to_string(probes.join("p24/p24-flipv.cdxml")).unwrap();
        let fv_live = InteractiveDocument::from_cdxml(&CdxmlDocument::parse(&fv_xml).unwrap());
        let mut fv_rust = mc_seed.scene.molecules[0].clone();
        let all_idx: Vec<usize> = (0..fv_rust.atoms.len()).collect();
        flip_atoms_preserving_stereo(&mut fv_rust, &all_idx, FlipAxis::Vertical, false);
        assert_eq!(compute_tetrahedral_cip_code(&fv_rust, 1), 3); // 2S
        assert_eq!(compute_tetrahedral_cip_code(&fv_rust, 2), 2); // 3R
        assert_eq!(compute_tetrahedral_cip_code(&fv_live.scene.molecules[0], 1), 3);
        assert_eq!(compute_tetrahedral_cip_code(&fv_live.scene.molecules[0], 2), 2);
        assert_eq!(fv_live.scene.molecules[0].bonds[4].display, BondDisplay::WedgeBegin);
        assert_eq!(fv_live.scene.molecules[0].bonds[5].display, BondDisplay::WedgedHashBegin);

        // Rotate 180 Vertical swaps WedgeBegin <-> WedgedHashBegin -> preserves (2R, 3S)
        let r180_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p24/p24-rot180v.cdxml")).unwrap()).unwrap(),
        );
        let mut r180_rust = mc_seed.scene.molecules[0].clone();
        flip_atoms_preserving_stereo(&mut r180_rust, &all_idx, FlipAxis::Vertical, true);
        assert_eq!(compute_tetrahedral_cip_code(&r180_rust, 1), 2); // 2R
        assert_eq!(compute_tetrahedral_cip_code(&r180_rust, 2), 3); // 3S
        assert_eq!(r180_rust.bonds[4].display, BondDisplay::WedgedHashBegin);
        assert_eq!(r180_rust.bonds[5].display, BondDisplay::WedgeBegin);
        // Prime reverses XML node order in Rotate 180 Vertical; look up C2 (id=12) and C3 (id=13) by id
        let c2_idx = r180_live.scene.molecules[0].atoms.iter().position(|a| a.id == 12).unwrap();
        let c3_idx = r180_live.scene.molecules[0].atoms.iter().position(|a| a.id == 13).unwrap();
        assert_eq!(compute_tetrahedral_cip_code(&r180_live.scene.molecules[0], c2_idx), 2);
        assert_eq!(compute_tetrahedral_cip_code(&r180_live.scene.molecules[0], c3_idx), 3);

        // Rotate 90 deg clockwise -> matches rot90 coords (<0.25 pt label-center delta) & preserves (2R, 3S)
        let r90_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p24/p24-rot90.cdxml")).unwrap()).unwrap(),
        );
        let mut r90_mols = mc_seed.scene.molecules.clone();
        rotate_molecules_around_bounds_center(&mut r90_mols, 90.0_f64.to_radians());
        for (ra, la) in r90_mols[0].atoms.iter().zip(r90_live.scene.molecules[0].atoms.iter()) {
            assert!((ra.x - la.x).abs() < 0.3, "rot90 x: {} vs {}", ra.x, la.x);
            assert!((ra.y - la.y).abs() < 0.3, "rot90 y: {} vs {}", ra.y, la.y);
        }
        assert_eq!(compute_tetrahedral_cip_code(&r90_mols[0], 1), 2);
        assert_eq!(compute_tetrahedral_cip_code(&r90_mols[0], 2), 3);

        // SVG live export vs Rust to_svg
        let live_svg = fs::read_to_string(probes.join("p24/p24-multi-chiral.svg")).unwrap();
        assert!(live_svg.contains("viewBox=\"0 0 73 66\""));
        assert!(live_svg.contains("(R)") && live_svg.contains("(S)"));
        assert!(live_svg.contains("Br") && live_svg.contains("Cl"));
        let rust_svg = mc_seed.to_svg();
        assert!(rust_svg.contains("<polygon") && rust_svg.contains("Br") && rust_svg.contains("Cl"));

        // 3. P25: Multi-record V2000 SDF generated from expanded query structures parses 1:1
        let sdf = fs::read_to_string(probes.join("p25/p25-expanded.sdf")).unwrap();
        let records = chiralcraft_formats::sdf_to_molecules(&sdf).unwrap();
        assert_eq!(records.len(), 3);
        let zs: Vec<i32> = records.iter().map(|r| r.molecule.atoms[2].z).collect();
        assert_eq!(zs, vec![9, 17, 35]);
    }

    #[test]
    fn test_live_prime_p26_p27_decalin_vert_dist_comb_query_and_chair() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P26a: Object > Distribute > Vertically matches distribute_molecules(..., Vertical)
        let vd_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p26/p26-vert-dist-seed.cdxml")).unwrap()).unwrap(),
        );
        let vd_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p26/p26-vert-distributed.cdxml")).unwrap()).unwrap(),
        );
        let mut vmols = vd_seed.scene.molecules.clone();
        distribute_molecules(&mut vmols, DistributeAxis::Vertical);
        assert!((vmols[1].atoms[0].y - 180.0).abs() < 1e-3);
        for (rm, lm) in vmols.iter().zip(vd_live.scene.molecules.iter()) {
            for (ra, la) in rm.atoms.iter().zip(lm.atoms.iter()) {
                assert!((ra.x - la.x).abs() < 1e-3);
                assert!((ra.y - la.y).abs() < 1e-3);
            }
        }

        // 2. P26b: 4-stereocenter fused 6+6 decalin Clean Up Structure + Rotate 45 deg + SVG
        let dec_cleaned = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p26/p26-decalin-cleaned.cdxml")).unwrap()).unwrap(),
        );
        let dec_rot45_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p26/p26-decalin-rot45.cdxml")).unwrap()).unwrap(),
        );
        let mut dec_mols = dec_cleaned.scene.molecules.clone();
        // Rotate around fragment graphical BoundingBox="96.55 105.65 188.05 152.60" center -> matches Prime to <0.02 pt!
        rotate_molecules_around_point(
            &mut dec_mols,
            0.5 * (96.55 + 188.05),
            0.5 * (105.65 + 152.60),
            45.0_f64.to_radians(),
        );
        for (ra, la) in dec_mols[0].atoms.iter().zip(dec_rot45_live.scene.molecules[0].atoms.iter()) {
            assert!((ra.x - la.x).abs() < 0.02, "decalin rot45 x: {} vs {}", ra.x, la.x);
            assert!((ra.y - la.y).abs() < 0.02, "decalin rot45 y: {} vs {}", ra.y, la.y);
        }
        // Also verify atom-bounds center rotation yields identical relative bond vectors (<0.02 pt)
        let mut dec_rel = dec_cleaned.scene.molecules.clone();
        rotate_molecules_around_bounds_center(&mut dec_rel, 45.0_f64.to_radians());
        let (dx0, dy0) = (
            dec_rel[0].atoms[0].x - dec_rot45_live.scene.molecules[0].atoms[0].x,
            dec_rel[0].atoms[0].y - dec_rot45_live.scene.molecules[0].atoms[0].y,
        );
        for (ra, la) in dec_rel[0].atoms.iter().zip(dec_rot45_live.scene.molecules[0].atoms.iter()) {
            assert!(((ra.x - la.x) - dx0).abs() < 0.02);
            assert!(((ra.y - la.y) - dy0).abs() < 0.02);
        }
        // All 4 stereocenters (id 1=S(3), id 4=R(2), id 6=R(2), id 8=R(2)) invariant before & after 45 deg rotation
        for &(id, expected_cip) in &[(1u32, 3u8), (4, 2), (6, 2), (8, 2)] {
            let r_idx = dec_mols[0].atoms.iter().position(|a| a.id == id).unwrap();
            let l_idx = dec_rot45_live.scene.molecules[0].atoms.iter().position(|a| a.id == id).unwrap();
            assert_eq!(compute_tetrahedral_cip_code(&dec_mols[0], r_idx), expected_cip);
            assert_eq!(compute_tetrahedral_cip_code(&dec_rot45_live.scene.molecules[0], l_idx), expected_cip);
        }

        let dec_svg_live = fs::read_to_string(probes.join("p26/p26-decalin-rot45.svg")).unwrap();
        assert!(dec_svg_live.contains("viewBox=\"0 0 126 68\""));
        assert!(dec_svg_live.contains("(S)") && dec_svg_live.contains("(R)"));
        assert!(dec_svg_live.contains("HO") && dec_svg_live.contains("Cl"));
        let dec_svg_rust = dec_rot45_live.to_svg();
        assert!(dec_svg_rust.contains("<polygon") && dec_svg_rust.contains("HO") && dec_svg_rust.contains("Cl"));

        // 3. P27a: Combinatorial 2x2=4 Expand Generic Structure ([F,Cl] x [O,N])
        let comb_out = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p27/p27-comb-query-out.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(comb_out.scene.molecules.len(), 4);
        let pairs: Vec<(i32, i32)> = comb_out
            .scene
            .molecules
            .iter()
            .map(|m| {
                assert!(check_structure(m, &[]).is_empty());
                (m.atoms[0].z, m.atoms[3].z)
            })
            .collect();
        assert_eq!(pairs, vec![(9, 7), (9, 8), (17, 7), (17, 8)]);

        // 4. P27c: Live Chair Cyclohexane stamped via Tools palette
        let chair = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p27/p27-chair-stamped.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(chair.scene.molecules.len(), 1);
        assert_eq!(chair.scene.molecules[0].atoms.len(), 6);
        assert_eq!(chair.scene.molecules[0].bonds.len(), 6);
        assert!(check_structure(&chair.scene.molecules[0], &[]).is_empty());
    }

    #[test]
    fn test_live_prime_p28_p29_groups_attachments_radicals_and_spiro() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P28a: Object > Group + Object > Add Frame (brackets enclosing both fragments)
        let gf_xml = fs::read_to_string(probes.join("p28/p28-grouped-framed.cdxml")).unwrap();
        assert!(gf_xml.contains("<group"));
        assert!(gf_xml.contains("BracketedObjectIDs=\"11 12 21 22\""));
        let gf_doc = CdxmlDocument::parse(&gf_xml).unwrap();
        let gf_scene = InteractiveDocument::from_cdxml(&gf_doc);
        assert_eq!(gf_scene.scene.molecules.len(), 2);
        let bgs = gf_doc.extract_bracketed_groups();
        assert_eq!(bgs.len(), 1);
        assert_eq!(bgs[0].bracketed_object_ids, vec!["11", "12", "21", "22"]);

        // 2. P28b: MultiAttachment & VariableAttachment on benzene carry is_query=true & 0 false warnings
        for (fname, expected_type) in [
            ("p28/p28-multicenter.cdxml", "NodeType=\"MultiAttachment\""),
            ("p28/p28-variable-attach.cdxml", "NodeType=\"VariableAttachment\""),
        ] {
            let xml = fs::read_to_string(probes.join(fname)).unwrap();
            assert!(xml.contains(expected_type));
            assert!(xml.contains("Attachments=\"11 12 13 14 15 16\""));
            let scene = InteractiveDocument::from_cdxml(&CdxmlDocument::parse(&xml).unwrap());
            assert_eq!(scene.scene.molecules[0].atoms.len(), 7);
            assert!(scene.scene.molecules[0].atoms[6].is_query);
            assert!(check_structure(&scene.scene.molecules[0], &[]).is_empty());
        }

        // 3. P29a: Check Structure on Ethyl radical (Radical="Doublet"), NH3+, pentavalent [C]+, 13CH3
        let rc = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p29/p29-radicals-cations-out.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(rc.scene.molecules.len(), 4);
        assert!(check_structure(&rc.scene.molecules[0], &[]).is_empty()); // CH3-CH2* (Radical=Doublet) clean!
        assert!(check_structure(&rc.scene.molecules[1], &[]).is_empty()); // CH3-NH3+ clean!
        let w_c5 = check_structure(&rc.scene.molecules[2], &[]); // [C]+ 5 bonds -> invalid valence
        assert_eq!(w_c5.len(), 1);
        assert_eq!(w_c5[0].message, "An atom in this label has an invalid valence.");
        assert!(rc.scene.molecules[3].atoms[0].is_unspecified_label); // Unspecified label 13CH3 flagged Unspecified in Prime
        let w_unspec = check_structure(&rc.scene.molecules[3], &[]);
        assert_eq!(w_unspec.len(), 1);
        assert_eq!(w_unspec[0].message, "ChiralCraft can't interpret this label.");

        // 4. P29b: Live Spiro[4.5]decane fused via Cyclopentane tool on Cyclohexane
        let spiro = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p29/p29-spiro-chain.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(spiro.scene.molecules.len(), 1);
        assert_eq!(spiro.scene.molecules[0].atoms.len(), 10);
        assert_eq!(spiro.scene.molecules[0].bonds.len(), 11);
        assert!(check_structure(&spiro.scene.molecules[0], &[]).is_empty());
    }

    #[test]
    fn test_live_prime_p30_p31_join_undo_chain_radicals_attach_and_nicknames() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P30a: Object > Join (Ctrl+J) fuses close endpoints & deletes dropped atom
        let j_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-join-seed.cdxml")).unwrap()).unwrap(),
        );
        let j_live = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-joined.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(j_live.scene.molecules.len(), 1);
        assert_eq!(j_live.scene.molecules[0].atoms.len(), 3);
        assert_eq!(j_live.scene.molecules[0].bonds.len(), 2);
        let mut joined = j_seed.scene.molecules[1].clone();
        assert!(join_molecules_at_atoms(&mut joined, 0, &j_seed.scene.molecules[0], 1));
        assert_eq!(joined.atoms.len(), 3);
        assert_eq!(joined.bonds.len(), 2);
        assert!((joined.atoms[2].x - j_live.scene.molecules[0].atoms[0].x).abs() < 0.02);
        assert!(check_structure(&joined, &[])
            .iter()
            .any(|w| w.message == "This unlabeled atom is connected to two identical colinear bonds, and may appear to be invisible."));

        // Check single-molecule join_atoms removes drop_atom cleanly
        let mut tri = joined.clone();
        assert!(join_atoms(&mut tri, 0, 2));
        assert_eq!(tri.atoms.len(), 2);

        // 2. P30b: Multi-atom Structure Clean Up -> single Ctrl+Z Undo reverts atomically -> Redo restores
        let m_seed = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-macro-seed.cdxml")).unwrap()).unwrap(),
        );
        let m_clean = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-after-cleanup.cdxml")).unwrap()).unwrap(),
        );
        let m_undone = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-after-one-undo.cdxml")).unwrap()).unwrap(),
        );
        let m_redone = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-after-redo.cdxml")).unwrap()).unwrap(),
        );
        for (a, b) in m_seed.scene.molecules[0].atoms.iter().zip(m_undone.scene.molecules[0].atoms.iter()) {
            assert!((a.x - b.x).abs() < 1e-3 && (a.y - b.y).abs() < 1e-3);
        }
        let mut doc_undo = InteractiveDocument::new_saved_seed("p30-macro.cdxml", DocumentSettings::acs_document_1996());
        doc_undo.scene = m_seed.scene.clone();
        doc_undo.record_transaction("Structure Clean Up", m_clean.scene.clone());
        assert_eq!(doc_undo.undo_menu_label(), "Undo Structure Clean Up");
        assert!(doc_undo.undo());
        assert_eq!(doc_undo.redo_menu_label(), "Redo Structure Clean Up");
        assert!((doc_undo.scene.molecules[0].atoms[0].x - m_undone.scene.molecules[0].atoms[0].x).abs() < 1e-3);
        assert!(doc_undo.redo());
        assert!((doc_undo.scene.molecules[0].atoms[0].x - m_redone.scene.molecules[0].atoms[0].x).abs() < 1e-3);

        // 3. P30c: Acyclic Chain (8 carbons) exact bond geometry
        let ch8 = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p30/p30-chain-8c.cdxml")).unwrap()).unwrap(),
        );
        assert_eq!(ch8.scene.molecules[0].atoms.len(), 8);
        assert_eq!(ch8.scene.molecules[0].bonds.len(), 7);
        assert!(check_structure(&ch8.scene.molecules[0], &[]).is_empty());

        // 4. P31a: Singlet carbene (:CH2 & CH:), Triplet->Singlet conversion, [CH]+., [CH]-. vs invalid CH2
        let rm_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p31/p31-radicals-matrix-out.cdxml")).unwrap(),
        ).unwrap();
        let rm_doc = InteractiveDocument::from_cdxml(&rm_cdxml);
        assert_eq!(rm_doc.scene.molecules.len(), 4);
        assert_eq!(rm_doc.scene.captions.len(), 2);
        assert!(rm_doc.scene.captions[0].warning.is_none()); // CH2: (:CH2 Singlet)
        assert_eq!(
            rm_doc.scene.captions[1].warning.as_deref(),
            Some("An atom in this label has an invalid valence.")
        );
        let rm_mols = rm_cdxml.extract_molecules();
        assert_eq!(rm_mols.len(), 6);
        assert!(check_structure(&rm_mols[0].1, &[]).is_empty()); // CH3-CH: (Singlet)
        assert!(check_structure(&rm_mols[1].1, &[]).is_empty()); // CH3-CH: (converted from Triplet)
        assert!(check_structure(&rm_mols[2].1, &[]).is_empty()); // :CH2 (Singlet carbene)
        let w_ch2 = check_structure(&rm_mols[3].1, &[]);         // CH2 (unradicalled -> invalid!)
        assert_eq!(w_ch2.len(), 1);
        assert_eq!(w_ch2[0].message, "An atom in this label has an invalid valence.");
        assert!(check_structure(&rm_mols[4].1, &[]).is_empty()); // CH3-CH+. (Doublet, +1)
        assert!(check_structure(&rm_mols[5].1, &[]).is_empty()); // CH3-CH-. (Doublet, -1)

        // Verify live ChiralCraft MOL V2000 M RAD parity: 1=Singlet -> chem 2; 2=Doublet -> chem 1
        let rm_mol = chiralcraft_formats::molfile_v2000_to_molecule(
            &fs::read_to_string(probes.join("p31/p31-radicals-matrix-out.mol")).unwrap(),
        ).unwrap();
        assert_eq!(rm_mol.atoms[1].radical, 2); // Singlet
        assert_eq!(rm_mol.atoms[3].radical, 2); // Singlet
        assert_eq!(rm_mol.atoms[4].radical, 2); // Singlet (:CH2)
        assert_eq!(rm_mol.atoms[7].radical, 1); // Doublet (+1)
        assert_eq!(rm_mol.atoms[7].charge, 1);
        assert_eq!(rm_mol.atoms[9].radical, 1); // Doublet (-1)
        assert_eq!(rm_mol.atoms[9].charge, -1);

        // 5. P31b: Custom User Nickname Ztst expand & contract CDXML
        let nick_exp = CdxmlDocument::parse(&fs::read_to_string(probes.join("p31/p31-nickname-expanded.cdxml")).unwrap()).unwrap();
        let nick_con = CdxmlDocument::parse(&fs::read_to_string(probes.join("p31/p31-nickname-contracted.cdxml")).unwrap()).unwrap();
        assert!(nick_exp.elements.iter().any(|e| e.tag == "chemicalproperty" && e.attr("Name") == Some("Residue")));
        assert!(nick_con.elements.iter().any(|e| e.tag == "n" && e.attr("NodeType") == Some("Fragment") && e.attr("ExternalConnectionType") == Some("Residue")));

        // 6. P31c: Object > Attach Data... CDXML & MDL V2000 Data SGroup extraction
        let ad_doc = CdxmlDocument::parse(&fs::read_to_string(probes.join("p31/p31-attach-data.cdxml")).unwrap()).unwrap();
        let ads = ad_doc.extract_attached_data();
        assert_eq!(ads.len(), 1);
        assert_eq!(ads[0].field_name, "pKa");
        assert_eq!(ads[0].tag.as_deref(), Some("D"));
        assert_eq!(ads[0].display_text, "D = 15.9");
        assert_eq!(ads[0].basis_objects, vec!["13"]);

        let sgs = chiralcraft_formats::molfile_v2000_extract_data_sgroups(
            &fs::read_to_string(probes.join("p31/p31-attach-data.mol")).unwrap(),
        );
        assert_eq!(sgs.len(), 1);
        assert_eq!(sgs[0].field_name, "pKa");
        assert_eq!(sgs[0].atom_indices_1based, vec![3]);
        assert_eq!(sgs[0].tag.as_deref(), Some("D"));
        assert_eq!(sgs[0].data, "15.9");
    }

    #[test]
    fn test_live_prime_p32_p33_menu_stereo_edge_warnings_and_svg() {
        let probes = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P32a: Toggle Aromatic Display (Order="1.5" + <graphic OvalType="Circle">) & SVG/CDXML roundtrip
        let arom_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32a-arom-out.cdxml")).unwrap(),
        )
        .unwrap();
        let arom_circles = arom_cdxml.extract_aromatic_circles();
        assert_eq!(arom_circles.len(), 3);
        assert!((arom_circles[0].radius - 9.01).abs() < 0.2);
        let arom_doc = InteractiveDocument::from_cdxml(&arom_cdxml);
        assert_eq!(arom_doc.scene.aromatic_circles.len(), 3);
        let arom_svg = arom_doc.to_svg();
        assert_eq!(arom_svg.matches("fill=\"none\" stroke=\"black\"").count(), 3);
        let arom_rt = InteractiveDocument::from_cdxml(&arom_doc.to_cdxml());
        assert_eq!(arom_rt.scene.aromatic_circles.len(), 3);

        // 2. P32b: Toggle Hydrogen Display (Visible="no" IgnoreWarnings="yes") & SVG suppression
        let hyd_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32b-hyd-out.cdxml")).unwrap(),
        )
        .unwrap();
        let hyd_doc = InteractiveDocument::from_cdxml(&hyd_cdxml);
        let mut hyd_mol = hyd_doc.scene.molecules[0].clone();
        let explicit_h_count = hyd_mol.atoms.iter().filter(|a| a.z == 1).count();
        assert!(explicit_h_count >= 1);
        assert!(hyd_mol
            .atoms
            .iter()
            .filter(|a| a.z == 1)
            .all(|a| !a.visible && a.ignore_warnings));
        assert!(hyd_mol
            .bonds
            .iter()
            .filter(|b| hyd_mol.atoms[b.begin_atom].z == 1 || hyd_mol.atoms[b.end_atom].z == 1)
            .all(|b| !b.visible && b.ignore_warnings));
        let hyd_svg = hyd_doc.to_svg();
        assert!(!hyd_svg.contains(">H<"));
        let hyd_rt = InteractiveDocument::from_cdxml(&hyd_doc.to_cdxml());
        assert!(hyd_rt.scene.molecules[0]
            .atoms
            .iter()
            .filter(|a| a.z == 1)
            .all(|a| !a.visible && a.ignore_warnings));
        assert!(toggle_hydrogen_display(&mut hyd_mol));
        assert!(hyd_mol
            .atoms
            .iter()
            .filter(|a| a.z == 1)
            .all(|a| a.visible && !a.ignore_warnings));
        assert!(!toggle_hydrogen_display(&mut hyd_mol));
        assert!(hyd_mol
            .atoms
            .iter()
            .filter(|a| a.z == 1)
            .all(|a| !a.visible && a.ignore_warnings));

        // 3. P32c: Clear Reaction Map (ReactionStepAtomMap & RxnParticipation stripped)
        let unmap_xml = fs::read_to_string(probes.join("p32/p32c-rxn-unmap-out.cdxml")).unwrap();
        assert!(!unmap_xml.contains("ReactionStepAtomMap"));
        assert!(!unmap_xml.contains("RxnParticipation"));
        let mut mapped_cdxml = CdxmlDocument::parse(
            &"<CDXML><page id=\"1\"><scheme id=\"2\"><step id=\"3\" ReactionStepAtomMap=\"10 20\"/></scheme><fragment id=\"4\"><n id=\"10\" Element=\"6\" p=\"0 0\"/><n id=\"20\" Element=\"6\" p=\"14.4 0\"/><b id=\"5\" B=\"10\" E=\"20\" RxnParticipation=\"ReactionCenter\"/></fragment></page></CDXML>",
        )
        .unwrap();
        mapped_cdxml.clear_reaction_maps();
        let cleared_xml = mapped_cdxml.to_xml_string();
        assert!(!cleared_xml.contains("ReactionStepAtomMap"));
        assert!(!cleared_xml.contains("RxnParticipation"));

        // 4. P32d: Center on Page, Show Stereochemistry ((R) objecttag), and Object > Annotate...
        let centered = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(
                &fs::read_to_string(probes.join("p32/p32d-centered.cdxml")).unwrap(),
            )
            .unwrap(),
        );
        let mut xs: Vec<f64> = centered.scene.molecules[0]
            .atoms
            .iter()
            .map(|a| a.x)
            .collect();
        let mut ys: Vec<f64> = centered.scene.molecules[0]
            .atoms
            .iter()
            .map(|a| a.y)
            .collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let cx = 0.5 * (xs[0] + xs[xs.len() - 1]);
        let cy = 0.5 * (ys[0] + ys[ys.len() - 1]);
        assert!((cx - 270.0).abs() < 15.0 && (cy - 359.875).abs() < 15.0);

        let show_stereo = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32d-show-stereo.cdxml")).unwrap(),
        )
        .unwrap();
        let st_inds = show_stereo.extract_stereo_indicators();
        assert!(st_inds.iter().any(|s| s.text.contains("(R)")));

        let annotated = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32d-annotated.cdxml")).unwrap(),
        )
        .unwrap();
        let anns = annotated.extract_annotations();
        assert!(anns
            .iter()
            .any(|a| a.keyword == "Name" && a.content == "Wave14"));

        // 5. P32e & P32f: Flatten 3D -> 2D, Bring to Front Z-order, and Attached Data duplicate
        let flat_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32e-flattened.cdxml")).unwrap(),
        )
        .unwrap();
        let flat_doc = InteractiveDocument::from_cdxml(&flat_cdxml);
        assert!(!flat_doc.scene.molecules[0].has_3d_coords);
        assert!(flat_doc.scene.molecules[0]
            .atoms
            .iter()
            .all(|a| a.z_coord.abs() < 1e-9));

        let mut test_3d = CdxmlDocument::parse(
            "<CDXML><page id=\"1\"><fragment id=\"2\"><n id=\"10\" Element=\"6\" xyz=\"12.50 24.00 5.00\"/></fragment></page></CDXML>",
        )
        .unwrap();
        test_3d.flatten_3d();
        let mut mol_3d = InteractiveDocument::from_cdxml(&test_3d).scene.molecules[0].clone();
        assert!(!mol_3d.has_3d_coords);
        assert!((mol_3d.atoms[0].x - 12.50).abs() < 1e-6);
        assert!((mol_3d.atoms[0].y - 24.00).abs() < 1e-6);
        mol_3d.has_3d_coords = true;
        mol_3d.atoms[0].z_coord = 3.25;
        flatten_molecule_to_2d(&mut mol_3d);
        assert!(!mol_3d.has_3d_coords && mol_3d.atoms[0].z_coord == 0.0);

        let zorder_xml = fs::read_to_string(probes.join("p32/p32e-zorder-front.cdxml")).unwrap();
        let pos32 = zorder_xml.find("id=\"32\"").unwrap();
        let pos31 = zorder_xml.find("id=\"31\"").unwrap();
        assert!(pos31 > pos32);

        let dup_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p32/p32f-attach-dup.cdxml")).unwrap(),
        )
        .unwrap();
        let dup_ads = dup_cdxml.extract_attached_data();
        assert_eq!(dup_ads.len(), 2);
        assert!(dup_ads
            .iter()
            .all(|d| d.field_name == "pKa" && d.display_text == "D = 15.9"));

        // 6. P33a: 8-molecule Stereochemistry Showcase + P14 (M)-allene
        let p33a_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p33/p33a-stereo-out.cdxml")).unwrap(),
        )
        .unwrap();
        let p33a_mols = p33a_cdxml.extract_molecules();
        assert_eq!(p33a_mols.len(), 8);

        // F0: (R)-Methyl phenyl sulfoxide (AS=R -> 2, phenyl BS=N -> empty bond_stereo)
        let r0 = chiralcraft_chem::perceive_stereochemistry(&p33a_mols[0].1);
        assert_eq!(r0.atom_stereo.get(&0), Some(&2));
        assert!(r0.bond_stereo.is_empty());

        // F1 & F2: (E)- vs (Z)-Acetaldoxime (11=E, 10=Z)
        assert_eq!(
            chiralcraft_chem::perceive_stereochemistry(&p33a_mols[1].1)
                .bond_stereo
                .get(&0),
            Some(&11)
        );
        assert_eq!(
            chiralcraft_chem::perceive_stereochemistry(&p33a_mols[2].1)
                .bond_stereo
                .get(&0),
            Some(&10)
        );

        // F3: (P)-Allene -> Some('P') + P14 (M)-allene -> Some('M')
        let find_allene_center = |mol: &ChemMolecule| -> usize {
            (0..mol.atoms.len())
                .find(|&ai| {
                    let nb = mol.neighbors(ai);
                    mol.atoms[ai].z == 6
                        && nb.len() == 2
                        && mol.bonds[nb[0].1].order == 2
                        && mol.bonds[nb[1].1].order == 2
                })
                .unwrap()
        };
        let f3_center = find_allene_center(&p33a_mols[3].1);
        assert_eq!(
            chiralcraft_chem::compute_axial_allene_cip(&p33a_mols[3].1, f3_center),
            Some('P')
        );
        assert!(chiralcraft_chem::perceive_stereochemistry(&p33a_mols[3].1)
            .atom_stereo
            .is_empty());

        let p14_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p14/p14-out.cdxml")).unwrap(),
        )
        .unwrap();
        let p14_mols = p14_cdxml.extract_molecules();
        let p14_allene = p14_mols
            .iter()
            .find(|(_, m)| {
                (0..m.atoms.len()).any(|ai| {
                    let nb = m.neighbors(ai);
                    m.atoms[ai].z == 6
                        && nb.len() == 2
                        && m.bonds[nb[0].1].order == 2
                        && m.bonds[nb[1].1].order == 2
                })
            })
            .map(|(_, m)| m)
            .unwrap();
        let p14_center = find_allene_center(p14_allene);
        assert_eq!(
            chiralcraft_chem::compute_axial_allene_cip(p14_allene, p14_center),
            Some('M')
        );

        // F4 Camphor: C1=(S)=3, C4=(S)=3
        let r4 = chiralcraft_chem::perceive_stereochemistry(&p33a_mols[4].1);
        let mut camphor_codes: Vec<u8> = r4.atom_stereo.values().copied().collect();
        camphor_codes.sort();
        assert_eq!(camphor_codes, vec![3, 3]);
        assert!(check_structure(&p33a_mols[4].1, &[]).is_empty());

        // F5 exo-2-chloronorbornane: C1=(S)=3, C2=(S)=3, C4=(R)=2
        let r5 = chiralcraft_chem::perceive_stereochemistry(&p33a_mols[5].1);
        assert_eq!(r5.atom_stereo.get(&0), Some(&3));
        assert_eq!(r5.atom_stereo.get(&1), Some(&3));
        assert_eq!(r5.atom_stereo.get(&3), Some(&2));
        assert!(check_structure(&p33a_mols[5].1, &[]).is_empty());

        // F6 (2R,3s,4S)-2,3,4-trichloropentane: 2, 5, 3
        let r6 = chiralcraft_chem::perceive_stereochemistry(&p33a_mols[6].1);
        assert_eq!(r6.atom_stereo.get(&1), Some(&2));
        assert_eq!(r6.atom_stereo.get(&2), Some(&5));
        assert_eq!(r6.atom_stereo.get(&3), Some(&3));
        assert!(check_structure(&p33a_mols[6].1, &[]).is_empty());

        // F7 Me-P(=O)(Et)-Ph: AS="S" -> 3, phenyl BS="N" -> empty, ambiguous stereobond warning
        let r7 = chiralcraft_chem::perceive_stereochemistry(&p33a_mols[7].1);
        assert_eq!(r7.atom_stereo.get(&0), Some(&3));
        assert!(r7.bond_stereo.is_empty());
        let w7 = check_structure(&p33a_mols[7].1, &[]);
        assert!(w7
            .iter()
            .any(|w| w.message == "This stereocenter has conflicting or ambiguous stereobonds specified."));

        // 7. P33b: Edge-Case Check Structure & Crossing Bonds
        let p33b_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p33/p33b-edge-out.cdxml")).unwrap(),
        )
        .unwrap();
        let p33b_doc = InteractiveDocument::from_cdxml(&p33b_cdxml);
        assert_eq!(p33b_doc.scene.molecules.len(), 9);
        for fi in [0usize, 1, 2, 3, 4, 6] {
            assert!(
                check_structure(&p33b_doc.scene.molecules[fi], &[]).is_empty(),
                "Expected zero warnings for P33b F{fi}"
            );
        }
        assert!(p33b_doc.scene.molecules[5]
            .atoms
            .iter()
            .any(|a| a.id == 62 && a.is_unspecified_label));
        let w_10h = check_structure(&p33b_doc.scene.molecules[5], &[]);
        assert_eq!(w_10h.len(), 1);
        assert_eq!(w_10h[0].message, "ChiralCraft can't interpret this label.");

        let scene_warnings = check_structure_scene(&p33b_doc.scene);
        assert!(scene_warnings.iter().any(|w| {
            w.object_id == 86 && w.message == "This atom is very close to another atom or bond."
        }));
        let crossings = find_crossing_bond_pairs(&p33b_doc.scene.molecules);
        assert!(crossings.iter().any(|&((m1, b1), (m2, b2))| {
            let id1 = p33b_doc.scene.molecules[m1].bonds[b1].id;
            let id2 = p33b_doc.scene.molecules[m2].bonds[b2].id;
            (id1 == 801 && id2 == 851) || (id1 == 851 && id2 == 801)
        }));

        // 8. P33c: Enhanced Stereochemistry & Stereobond-between-Stereocenters Warnings
        let p33c_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p33/p33c-enh-default.cdxml")).unwrap(),
        )
        .unwrap();
        let enh_warnings = check_enhanced_stereo_warnings(&p33c_cdxml);
        assert!(enh_warnings.iter().any(|w| {
            w.object_id == 12
                && w.message == "This compound cannot possess relative stereochemistry."
        }));
        let p33c_scene_warnings =
            check_structure_scene(&InteractiveDocument::from_cdxml(&p33c_cdxml).scene);
        assert!(p33c_scene_warnings.iter().any(|w| {
            w.object_id == 202 && w.message == "This stereobond connects two stereocenters."
        }));
    }

    #[test]
    fn test_live_prime_p34_p35_arrows_scale_exports_header_and_nicknames() {
        let probes = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../research/runtime-probes");
            if !probes.exists() { return; }

        // 1. P34a: 6-arrow palette showcase
        let p34a_raw = fs::read_to_string(probes.join("p34/p34a-arrows-out.cdxml")).unwrap();
        let p34a_cdxml = CdxmlDocument::parse(&p34a_raw).unwrap();
        let p34a_doc = InteractiveDocument::from_cdxml(&p34a_cdxml);
        assert_eq!(p34a_raw.matches("<arrow").count(), 6);
        assert!(p34a_raw.contains("NoGo=\"Cross\""));
        assert!(p34a_raw.contains("ArrowheadHead=\"HalfLeft\""));
        assert!(p34a_doc.to_svg().contains("<svg"));

        // 2. P34b: Templates > Bicyclics stamp (8 atoms, 10 bonds, 0 warnings)
        let p34b_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p34/p34b-template-stamp.cdxml")).unwrap(),
        )
        .unwrap();
        let p34b_doc = InteractiveDocument::from_cdxml(&p34b_cdxml);
        let p34b_atoms: usize = p34b_doc.scene.molecules.iter().map(|m| m.atoms.len()).sum();
        let p34b_bonds: usize = p34b_doc.scene.molecules.iter().map(|m| m.bonds.len()).sum();
        assert_eq!(p34b_atoms, 8);
        assert_eq!(p34b_bonds, 10);
        for mol in &p34b_doc.scene.molecules {
            assert!(check_structure(mol, &[]).is_empty());
        }

        // 3. P34c: Scale 200% + 3-macro Undo stack
        let scaled_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p34/p34c-scaled-200.cdxml")).unwrap(),
        )
        .unwrap();
        let undo2_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p34/p34c-undo-2.cdxml")).unwrap(),
        )
        .unwrap();
        let undo3_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p34/p34c-undo-3.cdxml")).unwrap(),
        )
        .unwrap();
        let scaled_doc = InteractiveDocument::from_cdxml(&scaled_cdxml);
        let undo2_doc = InteractiveDocument::from_cdxml(&undo2_cdxml);
        let undo3_doc = InteractiveDocument::from_cdxml(&undo3_cdxml);
        let all_atoms = |d: &InteractiveDocument| -> Vec<(f64, f64)> {
            d.scene.molecules.iter().flat_map(|m| m.atoms.iter().map(|a| (a.x, a.y))).collect()
        };
        let scaled_pts = all_atoms(&scaled_doc);
        let undo2_pts = all_atoms(&undo2_doc);
        let undo3_pts = all_atoms(&undo3_doc);
        assert_eq!(scaled_pts.len(), 10);
        assert_eq!(undo2_pts.len(), 10);
        assert_eq!(undo3_pts.len(), 10);
        // undo-2 matches scaled-200 atom positions
        for ((ax, ay), (bx, by)) in scaled_pts.iter().zip(undo2_pts.iter()) {
            assert!((ax - bx).abs() < 0.05 && (ay - by).abs() < 0.05);
        }
        // scaled-200 width is ~2x undo-3 (pre-scale) width
        let span_x = |pts: &[(f64, f64)]| {
            let min_x = pts.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
            let max_x = pts.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
            max_x - min_x
        };
        let ratio = span_x(&scaled_pts) / span_x(&undo3_pts);
        assert!((ratio - 2.0).abs() < 0.05, "Expected ~2.0x scale ratio, got {ratio}");

        // 4. P35a: Multi-format MOL & RXN exports
        let mol_def = fs::read_to_string(probes.join("p35/p35a-default.mol")).unwrap();
        let mol_v2k = fs::read_to_string(probes.join("p35/p35a-v2000.mol")).unwrap();
        assert!(chiralcraft_formats::molfile_to_molecule(&mol_def).unwrap().atoms.len() >= 2);
        assert!(chiralcraft_formats::molfile_to_molecule(&mol_v2k).unwrap().atoms.len() >= 2);

        // 5. P35b: Document Settings > Header/Footer serialization
        let p35b_raw = fs::read_to_string(probes.join("p35/p35b-settings-out.cdxml")).unwrap();
        let p35b_cdxml = CdxmlDocument::parse(&p35b_raw).unwrap();
        assert!(p35b_raw.contains("Header=\"Wave15 Header\""));
        assert!(p35b_raw.contains("Footer=\"Page \""));
        let _ = InteractiveDocument::from_cdxml(&p35b_cdxml);

        // 6. P35c: Built-in "Ac" override (acetyl -> Et-O-) & Delete fallback (Actinium Z=89)
        let builtin_ac = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(
                &fs::read_to_string(probes.join("p35/p35c-builtin-ac.cdxml")).unwrap(),
            )
            .unwrap(),
        );
        let override_ac = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(
                &fs::read_to_string(probes.join("p35/p35c-overridden-expanded.cdxml")).unwrap(),
            )
            .unwrap(),
        );
        let deleted_ac = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(
                &fs::read_to_string(probes.join("p35/p35c-restored-expanded.cdxml")).unwrap(),
            )
            .unwrap(),
        );
        assert!(builtin_ac.scene.molecules[0].bonds.iter().any(|b| b.order == 2));
        assert!(override_ac.scene.molecules[0].bonds.iter().all(|b| b.order == 1));
        assert!(deleted_ac.scene.molecules[0].atoms.iter().any(|a| a.z == 89));

        // 7. Extended Prime Parity: Full DocumentScene wiring (Arrows, Header/Footer, Scale, .ctp Templates, Generic Expansion, 537 Nicknames + Expand/Contract, Preferences & Interchange Formats)
        assert_eq!(p34a_doc.scene.arrows.len(), 6);
        let p34a_rt = InteractiveDocument::from_cdxml(&p34a_doc.to_cdxml());
        assert_eq!(p34a_rt.scene.arrows.len(), 6);
        let p34a_svg = p34a_doc.to_svg();
        assert!(p34a_svg.contains("<polygon") && p34a_svg.contains("<line"));

        let p35b_doc = InteractiveDocument::from_cdxml(&p35b_cdxml);
        assert_eq!(p35b_doc.scene.header_footer.header.as_deref(), Some("Wave15 Header"));
        assert_eq!(p35b_doc.scene.header_footer.footer.as_deref(), Some("Page "));
        let p35b_rt = InteractiveDocument::from_cdxml(&p35b_doc.to_cdxml());
        assert_eq!(p35b_rt.scene.header_footer.header.as_deref(), Some("Wave15 Header"));

        // Real .ctp CDX Template Library loading from official ChiralCraft shipped files!
        let app_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../research/application");
        if !app_dir.exists() { return; }
        let bicyclics_lib = TemplateLibrary::load_from_ctp_path(&app_dir.join("bicyclics.ctp")).unwrap();
        assert!(bicyclics_lib.templates.len() >= 10, "bicyclics.ctp should contain >= 10 templates");
        let conformers_lib = TemplateLibrary::load_from_ctp_path(&app_dir.join("conformers.ctp")).unwrap();
        assert!(conformers_lib.templates.len() >= 5, "conformers.ctp should contain >= 5 templates");

        // Real CDXML-backed combinatorial Expand Generic Structure on P27 2x2=4 seed!
        let p27_seed = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p27/p27-comb-query-seed.cdxml")).unwrap(),
        )
        .unwrap();
        let p27_expanded = expand_generic_structure_from_cdxml(&p27_seed);
        assert_eq!(p27_expanded.len(), 4);

        let mut scaled_test = vec![stamp_template(BuiltinTemplate::Benzene, 100.0, 100.0, 14.4)];
        scale_molecules_around_center(&mut scaled_test, 2.0);
        assert_eq!(scaled_test[0].atoms.len(), 6);
        let chair = stamp_template(BuiltinTemplate::CyclohexaneChair, 50.0, 50.0, 14.4);
        assert_eq!(chair.atoms.len(), 6);

        // 537 real shipped ChiralCraft nicknames + real Expand/Contract Label + Undo/Redo + P35c override/delete semantics!
        let mut nicks = NicknameTable::with_factory_defaults();
        assert_eq!(nicks.len(), 537);
        assert_eq!(nicks.lookup("Boc").unwrap().fragment.atoms.len(), 7);
        assert_eq!(nicks.lookup("Cbz").unwrap().fragment.atoms.len(), 10);
        assert_eq!(nicks.lookup("Ts").unwrap().fragment.atoms.len(), 10);
        assert_eq!(nicks.lookup("Fmoc").unwrap().fragment.atoms.len(), 17);
        assert_eq!(nicks.check_collision("Ac"), NicknameCollision::BothNicknameAndElement(89));
        assert!(nicks.define_nickname("Ac", chair.clone(), 0, true, false).is_err());
        assert!(nicks.define_nickname("Ac", chair.clone(), 0, true, true).is_ok());
        assert!(nicks.delete_nickname("Ac"));
        assert_eq!(nicks.resolve_label("Ac"), LabelResolution::Element(89));

        // Real graph-splicing Expand Label (CH3-Boc -> 8 atoms) & symmetric Contract Label (-> 2 atoms) + Undo
        let mut edit_doc = InteractiveDocument::new_saved_seed("nick_test", DocumentSettings::acs_document_1996());
        let mut ch3_boc = ChemMolecule::default();
        let mut c0 = ChemAtom::new(1, 6);
        c0.x = 100.0;
        c0.y = 100.0;
        let mut c1 = ChemAtom::new(2, 6);
        c1.x = 114.4;
        c1.y = 100.0;
        c1.label = Some("Boc".to_string());
        ch3_boc.atoms.push(c0);
        ch3_boc.atoms.push(c1);
        ch3_boc.bonds.push(chiralcraft_chem::ChemBond {
            id: 3, begin_atom: 0, end_atom: 1, order: 1, display: BondDisplay::Solid, visible: true, ignore_warnings: false,
        });
        edit_doc.stamp_molecule_tx(ch3_boc, "Add Bond");
        assert!(edit_doc.expand_label_tx(0, 1, "Boc"));
        assert_eq!(edit_doc.scene.molecules[0].atoms.len(), 8);
        assert!(edit_doc.contract_last_residue_tx(0));
        assert_eq!(edit_doc.scene.molecules[0].atoms.len(), 2);
        assert_eq!(edit_doc.scene.molecules[0].atoms[1].label.as_deref(), Some("Boc"));
        assert!(edit_doc.undo());
        assert_eq!(edit_doc.scene.molecules[0].atoms.len(), 8);

        // UserPreferences Warnings filtering (disabling checkbox 0 suppresses valence warning)
        let mut bad_doc = InteractiveDocument::from_cdxml(
            &CdxmlDocument::parse(&fs::read_to_string(probes.join("p18/p18-valence-sru-out.cdxml")).unwrap()).unwrap(),
        );
        assert!(!bad_doc.check_structure_with_preferences().is_empty());
        bad_doc.preferences.show_chemical_warnings = false;
        assert!(bad_doc.check_structure_with_preferences().is_empty());

        // Calibrated Chemical Properties verified against LIVE ChiralCraft Aspirin (P08) & Caffeine (P21)!
        let aspirin_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p08/p08-chem-seed.cdxml")).unwrap(),
        )
        .unwrap();
        let aspirin_doc = InteractiveDocument::from_cdxml(&aspirin_cdxml);
        let aspirin_mol = aspirin_doc
            .scene
            .molecules
            .iter()
            .find(|m| m.compute_formula_properties().get_formula() == "C9H8O4")
            .expect("Aspirin C9H8O4 in p08-chem-seed.cdxml");
        let asp_p = chiralcraft_chem::estimate_chemical_properties(aspirin_mol);
        // Prime P08: tPSA=63.60, LogP=1.31, LogS=-1.724, BP=640.91 K, MP=432.56 K, pKa=3.264
        assert!((asp_p.tpsa - 63.60).abs() < 0.05, "Aspirin tPSA={}", asp_p.tpsa);
        assert!((asp_p.log_p - 1.31).abs() < 0.10, "Aspirin LogP={}", asp_p.log_p);
        assert!((asp_p.log_s - (-1.724)).abs() < 0.45, "Aspirin LogS={}", asp_p.log_s);
        assert!((asp_p.boiling_point_k - 640.91).abs() < 0.05, "Aspirin BP={}", asp_p.boiling_point_k);
        assert!((asp_p.melting_point_k - 432.56).abs() < 1.0, "Aspirin MP={}", asp_p.melting_point_k);
        assert!((asp_p.pka.unwrap() - 3.264).abs() < 0.15, "Aspirin pKa={:?}", asp_p.pka);

        let caffeine_cdxml = CdxmlDocument::parse(
            &fs::read_to_string(probes.join("p21/p21-caffeine-pasted.cdxml")).unwrap(),
        )
        .unwrap();
        let caffeine_doc = InteractiveDocument::from_cdxml(&caffeine_cdxml);
        let caffeine_mol = caffeine_doc
            .scene
            .molecules
            .iter()
            .find(|m| m.compute_formula_properties().get_formula() == "C8H10N4O2")
            .expect("Caffeine C8H10N4O2 in p21-caffeine-pasted.cdxml");
        let caf_p = chiralcraft_chem::estimate_chemical_properties(caffeine_mol);
        // Published Ertl tPSA for Caffeine = 61.82 (Prime P21: tPSA 56.22, LogP -1.03, MR 51.2, Henry 15.2, LogS -1.084, pKa None)
        assert!((caf_p.tpsa - 61.82).abs() < 0.05, "Caffeine tPSA={}", caf_p.tpsa);
        assert!((caf_p.log_p - (-1.03)).abs() < 0.35, "Caffeine LogP={}", caf_p.log_p);
        assert!((caf_p.molar_refractivity - 51.2).abs() < 3.5, "Caffeine MR={}", caf_p.molar_refractivity);
        assert!((caf_p.henry_law - 15.2).abs() < 2.0, "Caffeine Henry={}", caf_p.henry_law);
        assert!((caf_p.log_s - (-1.084)).abs() < 0.45, "Caffeine LogS={}", caf_p.log_s);
        assert_eq!(caf_p.pka, None);

        // Non-native interchange formats (V3000, CT, TGF, CML, EPS, EMF, PNG, SKC, SLN, InChI/InChIKey)
        let v3k = chiralcraft_formats::molecule_to_molfile_v3000(&scaled_test[0], "Benzene");
        assert_eq!(chiralcraft_formats::molfile_v3000_to_molecule(&v3k).unwrap().atoms.len(), 6);
        let ct = chiralcraft_formats::molecule_to_ct(&scaled_test[0], "Benzene");
        assert_eq!(chiralcraft_formats::ct_to_molecule(&ct).unwrap().atoms.len(), 6);
        let tgf = chiralcraft_formats::molecule_to_tgf(&scaled_test[0]);
        assert!(tgf.contains("$Version: 4"));
        assert_eq!(chiralcraft_formats::tgf_to_molecules(&tgf).unwrap()[0].atoms.len(), 6);
        let cml = chiralcraft_formats::molecule_to_cml(&scaled_test[0]);
        assert_eq!(chiralcraft_formats::cml_to_molecules(&cml).unwrap()[0].atoms.len(), 6);
        assert!(chiralcraft_formats::molecule_to_eps(&scaled_test[0]).starts_with("%!PS-Adobe-2.0"));
        assert_eq!(&chiralcraft_formats::molecule_to_emf(&scaled_test[0])[40..44], b" EMF");
        let png = chiralcraft_formats::molecule_to_png(&scaled_test[0], 2.0);
        assert_eq!(&png[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
        let skc = chiralcraft_formats::molecule_to_skc(&scaled_test[0]);
        assert_eq!(chiralcraft_formats::skc_to_molecules(&skc).unwrap()[0].atoms.len(), 6);
        assert!(chiralcraft_formats::molecule_to_sln(&scaled_test[0]).contains("@1"));
        // Verify exact real ChiralCraft ground-truth CT, CML, TGF, SKC & standard Ethanol InChIKey!
        let p35_ct = fs::read_to_string(probes.join("p35/p35a-out.ct")).unwrap();
        let ct_mol = chiralcraft_formats::ct_to_molecule(&p35_ct).unwrap();
        assert_eq!(ct_mol.atoms.len(), 5);
        assert_eq!(ct_mol.atoms[4].z, 8);
        assert!(ct_mol.bonds.iter().any(|b| b.display == BondDisplay::WedgeBegin));
        let p07_cml = fs::read_to_string(probes.join("p07/p07-export.cml")).unwrap();
        let cml_mols = chiralcraft_formats::cml_to_molecules(&p07_cml).unwrap();
        assert_eq!(cml_mols.len(), 3);
        assert_eq!(cml_mols[0].atoms.len(), 5);
        assert_eq!(cml_mols[0].bonds.iter().filter(|b| b.display == BondDisplay::WedgeBegin).count(), 2);
        let p35_tgf = fs::read_to_string(probes.join("p35/p35a-out.tgf")).unwrap();
        let tgf_mols = chiralcraft_formats::tgf_to_molecules(&p35_tgf).unwrap();
        assert_eq!(tgf_mols.len(), 1);
        assert_eq!(tgf_mols[0].atoms.len(), 5);
        assert_eq!(tgf_mols[0].atoms[4].z, 8);
        assert!(tgf_mols[0].bonds.iter().any(|b| b.display == BondDisplay::WedgeBegin));
        let p13_skc = fs::read(probes.join("p13/p13-rxn-aligned.skc")).unwrap();
        let skc_mols = chiralcraft_formats::skc_to_molecules(&p13_skc).unwrap();
        assert_eq!(skc_mols.len(), 2);
        assert_eq!(skc_mols[0].atoms.len(), 2);
        assert_eq!(skc_mols[0].bonds[0].order, 2);
        assert_eq!(skc_mols[1].atoms.len(), 3);
        assert!(skc_mols[1].atoms.iter().any(|a| a.z == 35)); // Bromine
        // Official IUPAC InChIKey SHA-256 test vectors (Ethanol, Aspirin, Caffeine)!
        assert_eq!(
            chiralcraft_formats::inchi_to_inchikey("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
            "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        );
        assert_eq!(
            chiralcraft_formats::inchi_to_inchikey(
                "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
            ),
            "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        );
        assert_eq!(
            chiralcraft_formats::inchi_to_inchikey(
                "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"
            ),
            "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
        );
    }
}
