//! 1:1 Extended ChiralCraft Scene Objects, `.ctp` CDX Template Loader,
//! CDXML-backed Generic Structure Enumerator, Preferences Warning Filter,
//! CDXML Roundtrip Serialization, and Full Vector SVG Rendering (`P10`–`P35`).

use std::path::Path;

use chiralcraft_chem::{get_element_by_symbol, BondDisplay, ChemAtom, ChemBond, ChemMolecule};
use chiralcraft_formats::{
    CdxmlDocument, CdxmlNodeChild, CdxChild, CdxDocument, CdxObject, CDX_OBJ_GRAPHIC,
    CDX_OBJ_PAGE, CDX_PROP_BOUNDING_BOX,
};

#[derive(Debug, Clone, PartialEq)]
pub struct SceneArrow {
    pub id: u32,
    pub tail: (f64, f64),
    pub head: (f64, f64),
    pub arrowhead_head: String,
    pub arrowhead_tail: Option<String>,
    pub no_go: Option<String>,
    pub dipole: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SceneCurve {
    pub id: u32,
    pub points: Vec<(f64, f64)>,
    pub closed: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SceneGraphic {
    pub id: u32,
    pub graphic_type: String,
    pub sub_type: Option<String>,
    pub bbox: (f64, f64, f64, f64),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SceneTable {
    pub id: u32,
    pub rows: usize,
    pub cols: usize,
    pub bbox: (f64, f64, f64, f64),
    pub cells: Vec<Vec<String>>,
}

impl SceneTable {
    pub fn ensure_cells(&mut self) {
        self.rows = self.rows.max(1);
        self.cols = self.cols.max(1);
        self.cells.resize_with(self.rows, Vec::new);
        for row in &mut self.cells {
            row.resize_with(self.cols, String::new);
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SceneTlcPlate {
    pub id: u32,
    pub bbox: (f64, f64, f64, f64),
    pub origin_fraction: f64,
    pub spot_rfs: Vec<f64>,
    pub spot_labels: Vec<String>,
    pub spot_colors: Vec<(u8, u8, u8)>,
    pub solvent_label: String,
}

impl SceneTlcPlate {
    pub fn ensure_metadata(&mut self) {
        let n = self.spot_rfs.len();
        while self.spot_labels.len() < n {
            self.spot_labels.push(format!("L{}", self.spot_labels.len() + 1));
        }
        self.spot_labels.truncate(n);
        while self.spot_colors.len() < n {
            self.spot_colors.push((25, 25, 25));
        }
        self.spot_colors.truncate(n);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PageHeaderFooter {
    pub header: Option<String>,
    pub footer: Option<String>,
    pub width_pages: u32,
    pub height_pages: u32,
}

impl Default for PageHeaderFooter {
    fn default() -> Self {
        Self {
            header: None,
            footer: None,
            width_pages: 1,
            height_pages: 1,
        }
    }
}

/// ChiralCraft `File > Preferences...` (`Warnings` 13 toggles + `Building/Display` defaults, `P31`/`P33c`).
#[derive(Debug, Clone, PartialEq)]
pub struct UserPreferences {
    pub units: String,
    pub tolerance_px: u32,
    /// Factory defaults verified in `P33c`: `[true]*9 + [false] + [true]*3`.
    /// Indices:
    /// 0: Valence and Charge Errors
    /// 1: Unbalanced Parentheses
    /// 2: Invalid Isotopes
    /// 3: Isolated Bonds
    /// 4: Atoms Near Other Bonds
    /// 5: Unspecified Stereochemistry
    /// 6: Ambiguous Stereochemistry
    /// 7: Stereobond Between Chiral Atoms
    /// 8: Linear Atoms
    /// 9: Misc / Enhanced Stereo Relative Warnings (default OFF in Prime `P33c`)
    /// 10: Uninterpretable Labels
    /// 11: Crossing Bonds
    /// 12: Radical / Non-standard Valence
    pub warnings_enabled: [bool; 13],
    pub show_chemical_warnings: bool,
    pub fixed_bond_length_pt: f64,
}

impl Default for UserPreferences {
    fn default() -> Self {
        Self {
            units: "cm".to_string(),
            tolerance_px: 5,
            warnings_enabled: [
                true, true, true, true, true, true, true, true, true, false, true, true, true,
            ],
            show_chemical_warnings: true,
            fixed_bond_length_pt: 14.4,
        }
    }
}

impl UserPreferences {
    pub fn warning_index_for_message(msg: &str) -> Option<usize> {
        if msg.contains("invalid valence") || msg.contains("too many bonds") {
            Some(0)
        } else if msg.contains("Parentheses") {
            Some(1)
        } else if msg.contains("isolated bond") {
            Some(3)
        } else if msg.contains("very close to another atom or bond") {
            Some(4)
        } else if msg.contains("no stereobonds specified") {
            Some(5)
        } else if msg.contains("conflicting or ambiguous stereobonds") {
            Some(6)
        } else if msg.contains("stereobond connects two stereocenters") {
            Some(7)
        } else if msg.contains("colinear") || msg.contains("linear") {
            Some(8)
        } else if msg.contains("relative stereochemistry") {
            Some(9)
        } else if msg.contains("can't interpret this label") {
            Some(10)
        } else {
            None
        }
    }

    pub fn is_warning_allowed(&self, msg: &str) -> bool {
        if !self.show_chemical_warnings {
            return false;
        }
        match Self::warning_index_for_message(msg) {
            Some(idx) => self.warnings_enabled[idx],
            None => true,
        }
    }
}

/// Scale molecules around their combined bounding-box center (`Object > Scale...`, `P34c`).
pub fn scale_molecules_around_center(mols: &mut [ChemMolecule], factor: f64) {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for m in mols.iter() {
        for a in &m.atoms {
            min_x = min_x.min(a.x);
            min_y = min_y.min(a.y);
            max_x = max_x.max(a.x);
            max_y = max_y.max(a.y);
        }
    }
    if !min_x.is_finite() {
        return;
    }
    let cx = 0.5 * (min_x + max_x);
    let cy = 0.5 * (min_y + max_y);
    for m in mols.iter_mut() {
        for a in &mut m.atoms {
            a.x = cx + (a.x - cx) * factor;
            a.y = cy + (a.y - cy) * factor;
        }
    }
}

/// Combinatorial `Structure > Expand Generic Structure` (`P22`, `P25`, `P27`) directly from a parsed CDXML!
/// Reads `ElementList="..."` attributes on `<n>` elements (e.g. `ElementList="9 17"` or `"F Cl"` or `"8 7"`)
/// and enumerates the Cartesian product of concrete molecules.
pub fn expand_generic_structure_from_cdxml(cdxml: &CdxmlDocument) -> Vec<ChemMolecule> {
    let extracted = cdxml.extract_molecules();
    let Some((_occ, base_mol)) = extracted.into_iter().next() else {
        return Vec::new();
    };

    // Collect ElementList attributes keyed by node id
    let mut site_lists: Vec<(usize, Vec<i32>)> = Vec::new();
    for el in cdxml.elements.iter().filter(|e| e.tag == "n") {
        let Some(nid) = el.attr("id").and_then(|s| s.parse::<u32>().ok()) else {
            continue;
        };
        let Some(atom_idx) = base_mol.atoms.iter().position(|a| a.id == nid) else {
            continue;
        };
        if let Some(el_list) = el.attr("ElementList") {
            let tokens: Vec<&str> = el_list
                .split(|c: char| c.is_whitespace() || c == ',' || c == ';')
                .filter(|t| !t.is_empty())
                .collect();
            let is_not = tokens.first().map(|t| t.eq_ignore_ascii_case("NOT")).unwrap_or(false);
            let parsed: Vec<i32> = tokens
                .into_iter()
                .skip(usize::from(is_not))
                .filter_map(|tok| {
                    tok.parse::<i32>()
                        .ok()
                        .or_else(|| get_element_by_symbol(tok).map(|r| r.z))
                })
                .collect();
            let zs: Vec<i32> = if is_not {
                // Standard organic/medicinal complement set excluding negated elements
                [6, 7, 8, 9, 15, 16, 17, 35, 53]
                    .into_iter()
                    .filter(|z| !parsed.contains(z))
                    .collect()
            } else {
                parsed
            };
            if !zs.is_empty() {
                site_lists.push((atom_idx, zs));
            }
        }
    }
    expand_generic_structure(&base_mol, &site_lists)
}

/// Combinatorial `Structure > Expand Generic Structure` (`P22`, `P25`, `P27`), bounded to 512 products for safety.
pub fn expand_generic_structure(
    base: &ChemMolecule,
    site_element_lists: &[(usize, Vec<i32>)],
) -> Vec<ChemMolecule> {
    const MAX_EXPANDED: usize = 512;
    let mut results = vec![base.clone()];
    for (atom_idx, z_list) in site_element_lists {
        if z_list.is_empty() {
            continue;
        }
        let mut next = Vec::new();
        'outer: for mol in &results {
            for &z in z_list {
                let mut m = mol.clone();
                if let Some(a) = m.atoms.get_mut(*atom_idx) {
                    a.z = z;
                    a.is_query = false;
                }
                next.push(m);
                if next.len() >= MAX_EXPANDED {
                    break 'outer;
                }
            }
        }
        results = next;
    }
    results
}

/// Formats raw `.ctp`/`.cds` stem filenames into official ChiralCraft display names
/// (e.g., `clipware__part_1` -> `Clipware, Part 1`, `cp_rings` -> `CP Rings`, `dna_templates` -> `DNA Templates`).
pub fn pretty_template_asset_name(stem: &str) -> String {
    match stem {
        "advanced_biographics" => "Advanced BioGraphics".to_string(),
        "amino_acid_side_chains" => "Amino Acid Side Chains".to_string(),
        "amino_acids" => "Amino Acids".to_string(),
        "clipware__part_1" => "Clipware, Part 1".to_string(),
        "clipware__part_2" => "Clipware, Part 2".to_string(),
        "cp_rings" => "CP Rings".to_string(),
        "dna_templates" => "DNA Templates".to_string(),
        "ph_rings" => "Ph Rings".to_string(),
        "rna_templates" => "RNA Templates".to_string(),
        "acs_document_1996" => "ACS Document 1996".to_string(),
        "dark_mode" => "Dark Mode".to_string(),
        "i_draw_styles" => "iDraw Styles".to_string(),
        "j._het._chem" => "J. Het. Chem.".to_string(),
        "j._mol._mod.__1_column_" => "J. Mol. Mod. (1 column)".to_string(),
        "j._mol._mod.__2_column_" => "J. Mol. Mod. (2 column)".to_string(),
        "labeled_monomers" => "Labeled Monomers".to_string(),
        "nature_research_document" => "Nature Research Document".to_string(),
        "new_a4_document" => "New A4 Document".to_string(),
        "new_document" => "New Document".to_string(),
        "new_slide" => "New Slide".to_string(),
        "phytomedicine_document" => "Phytomedicine Document".to_string(),
        "rsc__1_column__document" => "RSC (1 column) Document".to_string(),
        "rsc__2_column__document" => "RSC (2 column) Document".to_string(),
        "science_of_synthesis" => "Science of Synthesis".to_string(),
        "wiley_document" => "Wiley Document".to_string(),
        other => other
            .split('_')
            .filter(|w| !w.is_empty())
            .map(|w| {
                let mut cs = w.chars();
                match cs.next() {
                    Some(c) => format!("{}{}", c.to_uppercase(), cs.as_str()),
                    None => String::new(),
                }
            })
            .collect::<Vec<_>>()
            .join(" "),
    }
}

/// Real `.ctp` CDX Template Library loader — parses ALL 26 shipped ChiralCraft `.ctp` CDX files,
/// loading both chemical fragment templates (`templates`) and vector graphic/curve templates (`graphic_templates`)!
#[derive(Debug, Clone)]
pub struct TemplateLibrary {
    pub name: String,
    pub templates: Vec<ChemMolecule>,
    pub graphic_templates: Vec<Vec<SceneCurve>>,
}

impl TemplateLibrary {
    pub fn load_from_ctp_bytes(name: &str, bytes: &[u8]) -> Result<Self, String> {
        let cdx = CdxDocument::parse(bytes).map_err(|e| format!("{e:?}"))?;
        let cdxml = cdx.to_cdxml();
        let templates: Vec<ChemMolecule> = cdxml
            .extract_molecules()
            .into_iter()
            .map(|(_, m)| m)
            .filter(|m| !m.atoms.is_empty())
            .collect();

        // Also extract per-page vector curve & graphic templates directly from the CDX tree
        // (`0x8008` kCDXObj_Curve + `0x0A23` kCDXProp_Curve_Points + `0x8007` kCDXObj_Graphic)!
        let mut graphic_templates: Vec<Vec<SceneCurve>> = Vec::new();
        let mut next_cid = 1_u32;

        fn collect_cdx_curves(obj: &CdxObject, out: &mut Vec<SceneCurve>, next_cid: &mut u32) {
            const CDX_OBJ_CURVE: u16 = 0x8008;
            const CDX_PROP_CURVE_TYPE: u16 = 0x0A08;
            const CDX_PROP_CURVE_POINTS: u16 = 0x0A23;
            const CDX_PROP_GRAPHIC_TYPE: u16 = 0x0A00;

            if obj.tag == CDX_OBJ_CURVE {
                let mut closed = false;
                let mut pts: Vec<(f64, f64)> = Vec::new();
                for ch in &obj.children {
                    if let CdxChild::Property(p) = ch {
                        if p.tag == CDX_PROP_CURVE_TYPE && p.payload.len() >= 2 {
                            let ctype = u16::from_le_bytes([p.payload[0], p.payload[1]]);
                            closed = (ctype & 0x0001) != 0;
                        } else if p.tag == CDX_PROP_CURVE_POINTS && p.payload.len() >= 2 {
                            let cnt = u16::from_le_bytes([p.payload[0], p.payload[1]]) as usize;
                            let avail = (p.payload.len() - 2) / 8;
                            let n = cnt.min(avail);
                            // Subsample cubic control handles (keep anchors & midpoints smoothly)
                            for i in 0..n {
                                let off = 2 + i * 8;
                                let y = i32::from_le_bytes([
                                    p.payload[off], p.payload[off + 1], p.payload[off + 2], p.payload[off + 3],
                                ]) as f64
                                    / 65536.0;
                                let x = i32::from_le_bytes([
                                    p.payload[off + 4], p.payload[off + 5], p.payload[off + 6], p.payload[off + 7],
                                ]) as f64
                                    / 65536.0;
                                if pts
                                    .last()
                                    .map(|&(px, py)| (px - x).hypot(py - y) > 0.25)
                                    .unwrap_or(true)
                                {
                                    pts.push((x, y));
                                }
                            }
                        }
                    }
                }
                if pts.len() >= 2 {
                    out.push(SceneCurve {
                        id: *next_cid,
                        points: pts,
                        closed,
                    });
                    *next_cid += 1;
                }
            } else if obj.tag == CDX_OBJ_GRAPHIC {
                let mut gtype = 1_u16;
                let mut bbox: Option<(f64, f64, f64, f64)> = None;
                for ch in &obj.children {
                    if let CdxChild::Property(p) = ch {
                        if p.tag == CDX_PROP_GRAPHIC_TYPE && p.payload.len() >= 2 {
                            gtype = u16::from_le_bytes([p.payload[0], p.payload[1]]);
                        } else if p.tag == CDX_PROP_BOUNDING_BOX && p.payload.len() >= 16 {
                            let top = i32::from_le_bytes([p.payload[0], p.payload[1], p.payload[2], p.payload[3]]) as f64 / 65536.0;
                            let left = i32::from_le_bytes([p.payload[4], p.payload[5], p.payload[6], p.payload[7]]) as f64 / 65536.0;
                            let bottom = i32::from_le_bytes([p.payload[8], p.payload[9], p.payload[10], p.payload[11]]) as f64 / 65536.0;
                            let right = i32::from_le_bytes([p.payload[12], p.payload[13], p.payload[14], p.payload[15]]) as f64 / 65536.0;
                            bbox = Some((left, top, right, bottom));
                        }
                    }
                }
                if let Some((x1, y1, x2, y2)) = bbox {
                    let (pts, closed) = match gtype {
                        3 => (vec![(x1, y1), (x2, y1), (x2, y2), (x1, y2)], true),
                        4 => {
                            let cx = 0.5 * (x1 + x2);
                            let cy = 0.5 * (y1 + y2);
                            let rx = 0.5 * (x2 - x1).abs();
                            let ry = 0.5 * (y2 - y1).abs();
                            let ring = (0..12)
                                .map(|k| {
                                    let a = (k as f64) * std::f64::consts::TAU / 12.0;
                                    (cx + rx * a.cos(), cy + ry * a.sin())
                                })
                                .collect();
                            (ring, true)
                        }
                        _ => (vec![(x1, y1), (x2, y2)], false),
                    };
                    if pts.len() >= 2 {
                        out.push(SceneCurve {
                            id: *next_cid,
                            points: pts,
                            closed,
                        });
                        *next_cid += 1;
                    }
                }
            }
            for ch in &obj.children {
                if let CdxChild::Object(sub) = ch {
                    collect_cdx_curves(sub, out, next_cid);
                }
            }
        }

        if let Some(root) = cdx.root_objects.first() {
            for ch in &root.children {
                if let CdxChild::Object(page) = ch {
                    if page.tag == CDX_OBJ_PAGE {
                        let mut page_curves = Vec::new();
                        collect_cdx_curves(page, &mut page_curves, &mut next_cid);
                        if !page_curves.is_empty() {
                            graphic_templates.push(page_curves);
                        }
                    }
                }
            }
        }

        Ok(Self {
            name: pretty_template_asset_name(name),
            templates,
            graphic_templates,
        })
    }

    pub fn total_count(&self) -> usize {
        if !self.templates.is_empty() {
            self.templates.len()
        } else {
            self.graphic_templates.len()
        }
    }

    /// 100% clean-room built-in chemical & bio-graphics template libraries generated purely in Rust
    /// (zero external or proprietary `.ctp` assets required).
    pub fn builtin_libraries() -> Vec<Self> {
        let from_smiles_list = |name: &str, list: &[&str]| -> Self {
            let mut templates = Vec::new();
            for &smi in list {
                if let Ok(mut m) = chiralcraft_formats::smiles_to_molecule(smi) {
                    crate::clean_up_structure_preserving_stereo(&mut m, &[], 14.40);
                    if !m.atoms.is_empty() {
                        templates.push(m);
                    }
                }
            }
            Self {
                name: name.to_string(),
                templates,
                graphic_templates: Vec::new(),
            }
        };

        let mut libs = vec![
            from_smiles_list(
                "Aromatics & Heterocycles",
                &[
                    "c1ccccc1",
                    "c1ccc2ccccc2c1",
                    "c1ccc2cc3ccccc3cc2c1",
                    "c1ccc2c(c1)ccc1ccccc12",
                    "c1ccc(-c2ccccc2)cc1",
                    "c1ccncc1",
                    "c1cncnc1",
                    "c1cnccn1",
                    "c1cc[nH]c1",
                    "c1ccoc1",
                    "c1ccsc1",
                    "c1c[nH]cn1",
                    "c1cocn1",
                    "c1cscn1",
                    "c1ccc2[nH]ccc2c1",
                    "c1ccc2ncccc2c1",
                    "c1ccc2cnccc2c1",
                    "c1ncc2[nH]cnc2n1",
                ],
            ),
            from_smiles_list(
                "Amino Acids",
                &[
                    "NCC(=O)O",
                    "CC(N)C(=O)O",
                    "CC(C)C(N)C(=O)O",
                    "CC(C)CC(N)C(=O)O",
                    "CCC(C)C(N)C(=O)O",
                    "NC(CO)C(=O)O",
                    "CC(O)C(N)C(=O)O",
                    "NC(CS)C(=O)O",
                    "CSCCC(N)C(=O)O",
                    "NC(CC(=O)O)C(=O)O",
                    "NC(CCC(=O)O)C(=O)O",
                    "NC(CC(=O)N)C(=O)O",
                    "NC(CCC(=O)N)C(=O)O",
                    "NCCCC(N)C(=O)O",
                    "NC(CCCNC(=N)N)C(=O)O",
                    "NC(Cc1c[nH]cn1)C(=O)O",
                    "NC(Cc1ccccc1)C(=O)O",
                    "NC(Cc1ccc(O)cc1)C(=O)O",
                    "NC(Cc1c[nH]c2ccccc12)C(=O)O",
                    "O=C(O)C1CCCN1",
                ],
            ),
            from_smiles_list(
                "Rings & Polycycles",
                &[
                    "C1CC1",
                    "C1CCC1",
                    "C1CCCC1",
                    "C1CCCCC1",
                    "C1CCCCCC1",
                    "C1CCCCCCC1",
                    "C1=CCC=C1",
                    "C1CCNCC1",
                    "C1CNCCN1",
                    "C1COCCN1",
                    "C1CCNC1",
                    "C1CCOC1",
                    "C1COCCO1",
                    "C1CCC2CCCCC2C1",
                    "C1CC2CCC1C2",
                    "C1C2CC3CC1CC(C2)C3",
                ],
            ),
            from_smiles_list(
                "Carbohydrates",
                &[
                    "OCC1OC(O)C(O)C(O)C1O",
                    "OCC1OC(O)C(O)C1O",
                    "OCC1OC(O)CC1O",
                    "OCC1(O)OC(CO)C(O)C1O",
                    "OCC(O)C(O)C(O)C(O)C=O",
                    "C1CCOC1",
                    "C1CCOCC1",
                ],
            ),
            from_smiles_list(
                "Nucleobases",
                &[
                    "Nc1ncnc2[nH]cnc12",
                    "Nc1nc2[nH]cnc2c(=O)[nH]1",
                    "Nc1cc[nH]c(=O)n1",
                    "Cc1c[nH]c(=O)[nH]c1=O",
                    "O=c1cc[nH]c(=O)[nH]1",
                ],
            ),
            from_smiles_list(
                "Pharmaceuticals",
                &[
                    "CC(=O)Oc1ccccc1C(=O)O",
                    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
                    "CC(=O)Nc1ccc(O)cc1",
                    "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
                    "NCCc1ccc(O)c(O)c1",
                    "NCCc1c[nH]c2ccc(O)cc12",
                    "CN1CCCC1c1cccnc1",
                    "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
                ],
            ),
            from_smiles_list(
                "Functional & Protecting Groups",
                &[
                    "CC(=O)O",
                    "CCOC(=O)C",
                    "CC(=O)N",
                    "CS(=O)(=O)N",
                    "CC(C)(C)OC(=O)N",
                    "Cc1ccc(S(=O)(=O)O)cc1",
                    "FC(F)(F)S(=O)(=O)O",
                    "O=[N+]([O-])c1ccccc1",
                    "CC#N",
                    "OB(O)c1ccccc1",
                    "COP(=O)(OC)OC",
                ],
            ),
        ];

        // Clean-room vector BioGraphics & Lab Shapes (pure geometry, zero external files)
        let circle = |id: u32, cx: f64, cy: f64, rx: f64, ry: f64, n: usize| -> SceneCurve {
            let pts = (0..n)
                .map(|k| {
                    let a = (k as f64) * std::f64::consts::TAU / (n as f64);
                    (cx + rx * a.cos(), cy + ry * a.sin())
                })
                .collect();
            SceneCurve { id, points: pts, closed: true }
        };
        let mut dna_a = Vec::new();
        let mut dna_b = Vec::new();
        for k in 0..=18 {
            let t = (k as f64) / 18.0;
            let x = t * 64.0;
            let s = (t * std::f64::consts::TAU * 1.5).sin() * 10.0;
            dna_a.push((x, 20.0 + s));
            dna_b.push((x, 20.0 - s));
        }
        let bio_graphics: Vec<Vec<SceneCurve>> = vec![
            // 1. Cell / Vesicle
            vec![circle(1, 0.0, 0.0, 26.0, 20.0, 20), circle(2, 4.0, -2.0, 9.0, 7.0, 14)],
            // 2. Mitochondrion
            vec![
                circle(1, 0.0, 0.0, 30.0, 14.0, 20),
                SceneCurve {
                    id: 2,
                    points: vec![(-20.0, 0.0), (-14.0, -7.0), (-8.0, 7.0), (-2.0, -7.0), (4.0, 7.0), (10.0, -7.0), (16.0, 7.0), (20.0, 0.0)],
                    closed: false,
                },
            ],
            // 3. Y-shaped Antibody / Receptor
            vec![
                SceneCurve { id: 1, points: vec![(0.0, 24.0), (0.0, 0.0), (-16.0, -18.0)], closed: false },
                SceneCurve { id: 2, points: vec![(0.0, 0.0), (16.0, -18.0)], closed: false },
                SceneCurve { id: 3, points: vec![(-5.0, -6.0), (-13.0, -15.0)], closed: false },
                SceneCurve { id: 4, points: vec![(5.0, -6.0), (13.0, -15.0)], closed: false },
            ],
            // 4. DNA Double Helix
            vec![
                SceneCurve { id: 1, points: dna_a, closed: false },
                SceneCurve { id: 2, points: dna_b, closed: false },
            ],
            // 5. Lipid Bilayer Segment
            vec![
                SceneCurve { id: 1, points: vec![(-30.0, -8.0), (30.0, -8.0)], closed: false },
                SceneCurve { id: 2, points: vec![(-30.0, 8.0), (30.0, 8.0)], closed: false },
                circle(3, -20.0, -12.0, 4.0, 4.0, 10),
                circle(4, 0.0, -12.0, 4.0, 4.0, 10),
                circle(5, 20.0, -12.0, 4.0, 4.0, 10),
                circle(6, -20.0, 12.0, 4.0, 4.0, 10),
                circle(7, 0.0, 12.0, 4.0, 4.0, 10),
                circle(8, 20.0, 12.0, 4.0, 4.0, 10),
            ],
            // 6. Erlenmeyer Flask
            vec![SceneCurve {
                id: 1,
                points: vec![(-7.0, -24.0), (7.0, -24.0), (6.0, -8.0), (22.0, 22.0), (-22.0, 22.0), (-6.0, -8.0)],
                closed: true,
            }],
            // 7. Round-Bottom Flask
            vec![
                circle(1, 0.0, 8.0, 18.0, 18.0, 20),
                SceneCurve { id: 2, points: vec![(-6.0, -24.0), (-6.0, -9.0), (6.0, -9.0), (6.0, -24.0)], closed: false },
            ],
            // 8. Enzyme / Active-Site Pacman
            vec![SceneCurve {
                id: 1,
                points: {
                    let mut pts = vec![(0.0, 0.0)];
                    for k in 2..=14 {
                        let a = (k as f64) * std::f64::consts::TAU / 16.0;
                        pts.push((22.0 * a.cos(), 22.0 * a.sin()));
                    }
                    pts
                },
                closed: true,
            }],
        ];
        libs.push(Self {
            name: "BioGraphics & Lab Shapes".to_string(),
            templates: Vec::new(),
            graphic_templates: bio_graphics,
        });

        // Library 9: Dedicated pure-geometry Haworth Sugar Projections, Macrocycles & Bridged Cages!
        libs.push(Self {
            name: "Haworth & Macrocycles".to_string(),
            templates: vec![
                build_haworth_pyranose(14.40),
                build_haworth_furanose(14.40),
                build_porphine_core(14.40),
                build_crown_ether_18_6(14.40),
                build_adamantane_2d(14.40),
                build_norbornane_2d(14.40),
            ],
            graphic_templates: Vec::new(),
        });

        // Library 10: Modular Lab Glassware & Gas-Line Apparatus (pure vector schematics)
        let glassware: Vec<Vec<SceneCurve>> = vec![
            // 1. 1-Neck Round-Bottom Flask
            vec![
                circle(1, 0.0, 12.0, 20.0, 20.0, 24),
                SceneCurve { id: 2, points: vec![(-6.0, -26.0), (-6.0, -7.0), (6.0, -7.0), (6.0, -26.0)], closed: false },
                SceneCurve { id: 3, points: vec![(-8.0, -26.0), (8.0, -26.0)], closed: false },
            ],
            // 2. 3-Neck Round-Bottom Flask
            vec![
                circle(1, 0.0, 12.0, 22.0, 20.0, 24),
                SceneCurve { id: 2, points: vec![(-5.0, -26.0), (-5.0, -8.0), (5.0, -8.0), (5.0, -26.0)], closed: false },
                SceneCurve { id: 3, points: vec![(-22.0, -14.0), (-14.0, -2.0), (-8.0, -6.0), (-16.0, -18.0)], closed: false },
                SceneCurve { id: 4, points: vec![(22.0, -14.0), (14.0, -2.0), (8.0, -6.0), (16.0, -18.0)], closed: false },
            ],
            // 3. Liebig Reflux Condenser
            vec![
                SceneCurve { id: 1, points: vec![(-4.0, -32.0), (-4.0, 32.0)], closed: false },
                SceneCurve { id: 2, points: vec![(4.0, -32.0), (4.0, 32.0)], closed: false },
                SceneCurve { id: 3, points: vec![(-11.0, -20.0), (11.0, -20.0), (11.0, 20.0), (-11.0, 20.0)], closed: true },
                SceneCurve { id: 4, points: vec![(-11.0, 14.0), (-20.0, 14.0)], closed: false },
                SceneCurve { id: 5, points: vec![(11.0, -14.0), (20.0, -14.0)], closed: false },
            ],
            // 4. Separatory / Dropping Funnel
            vec![
                SceneCurve { id: 1, points: vec![(-5.0, -28.0), (5.0, -28.0), (18.0, -8.0), (3.0, 18.0), (3.0, 30.0), (-3.0, 30.0), (-3.0, 18.0), (-18.0, -8.0)], closed: true },
                SceneCurve { id: 2, points: vec![(-9.0, 21.0), (9.0, 21.0)], closed: false },
            ],
            // 5. Schlenk Vacuum / Gas Manifold (3-Port)
            vec![
                SceneCurve { id: 1, points: vec![(-36.0, -10.0), (36.0, -10.0), (36.0, -4.0), (-36.0, -4.0)], closed: true },
                SceneCurve { id: 2, points: vec![(-36.0, 4.0), (36.0, 4.0), (36.0, 10.0), (-36.0, 10.0)], closed: true },
                SceneCurve { id: 3, points: vec![(-22.0, -4.0), (-22.0, 22.0)], closed: false },
                SceneCurve { id: 4, points: vec![(0.0, -4.0), (0.0, 22.0)], closed: false },
                SceneCurve { id: 5, points: vec![(22.0, -4.0), (22.0, 22.0)], closed: false },
                circle(6, -22.0, 14.0, 4.0, 4.0, 10),
                circle(7, 0.0, 14.0, 4.0, 4.0, 10),
                circle(8, 22.0, 14.0, 4.0, 4.0, 10),
            ],
            // 6. Tube Furnace
            vec![
                SceneCurve { id: 1, points: vec![(-24.0, -14.0), (24.0, -14.0), (24.0, 14.0), (-24.0, 14.0)], closed: true },
                SceneCurve { id: 2, points: vec![(-38.0, -4.0), (38.0, -4.0)], closed: false },
                SceneCurve { id: 3, points: vec![(-38.0, 4.0), (38.0, 4.0)], closed: false },
                SceneCurve { id: 4, points: vec![(-16.0, -10.0), (-8.0, -6.0), (0.0, -10.0), (8.0, -6.0), (16.0, -10.0)], closed: false },
            ],
            // 7. Mass Flow Controller (MFC)
            vec![
                SceneCurve { id: 1, points: vec![(-18.0, -14.0), (18.0, -14.0), (18.0, 14.0), (-18.0, 14.0)], closed: true },
                SceneCurve { id: 2, points: vec![(-30.0, 0.0), (-18.0, 0.0)], closed: false },
                SceneCurve { id: 3, points: vec![(18.0, 0.0), (30.0, 0.0)], closed: false },
                SceneCurve { id: 4, points: vec![(-10.0, 6.0), (0.0, -6.0), (10.0, 6.0)], closed: false },
            ],
            // 8. 3-Way Valve / Gas-Line Tee + Elbow Tubing
            vec![
                SceneCurve { id: 1, points: vec![(-24.0, 0.0), (-8.0, -8.0), (-8.0, 8.0)], closed: true },
                SceneCurve { id: 2, points: vec![(24.0, 0.0), (8.0, -8.0), (8.0, 8.0)], closed: true },
                SceneCurve { id: 3, points: vec![(0.0, 22.0), (-8.0, 8.0), (8.0, 8.0)], closed: true },
                circle(4, 0.0, 0.0, 5.0, 5.0, 12),
            ],
        ];
        libs.push(Self {
            name: "Lab Glassware & Apparatus".to_string(),
            templates: Vec::new(),
            graphic_templates: glassware,
        });
        // Library 11: Clean-Room Organometallic Catalysts
        libs.push(Self {
            name: "Organometallic Catalysts".to_string(),
            templates: vec![
                crate::build_pd_pph3_4(14.40),
                crate::build_cisplatin(14.40),
                crate::build_wilkinson_catalyst(14.40),
                crate::build_vaska_complex(14.40),
                crate::build_ferrocene_sandwich(14.40),
                crate::build_grubbs_core(14.40),
            ],
            graphic_templates: Vec::new(),
        });
        libs
    }

    pub fn stamp_graphic(&self, idx: usize, cx: f64, cy: f64) -> Option<Vec<SceneCurve>> {
        let mut curves = self.graphic_templates.get(idx)?.clone();
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for c in &curves {
            for &(x, y) in &c.points {
                min_x = min_x.min(x);
                min_y = min_y.min(y);
                max_x = max_x.max(x);
                max_y = max_y.max(y);
            }
        }
        if min_x.is_finite() {
            let w = (max_x - min_x).max(1.0);
            let h = (max_y - min_y).max(1.0);
            let scale = (72.0 / w.max(h)).min(1.0);
            let ox = 0.5 * (min_x + max_x);
            let oy = 0.5 * (min_y + max_y);
            for c in &mut curves {
                for p in &mut c.points {
                    p.0 = cx + (p.0 - ox) * scale;
                    p.1 = cy + (p.1 - oy) * scale;
                }
            }
        }
        Some(curves)
    }

    pub fn load_from_ctp_path(path: &Path) -> Result<Self, String> {
        let bytes = std::fs::read(path).map_err(|e| e.to_string())?;
        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("templates");
        Self::load_from_ctp_bytes(name, &bytes)
    }

    /// Translate template `idx` so its bounding-box center lands at `(cx, cy)`.
    pub fn stamp(&self, idx: usize, cx: f64, cy: f64) -> Option<ChemMolecule> {
        let mut m = self.templates.get(idx)?.clone();
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for a in &m.atoms {
            min_x = min_x.min(a.x);
            min_y = min_y.min(a.y);
            max_x = max_x.max(a.x);
            max_y = max_y.max(a.y);
        }
        if min_x.is_finite() {
            let dx = cx - 0.5 * (min_x + max_x);
            let dy = cy - 0.5 * (min_y + max_y);
            for a in &mut m.atoms {
                a.x += dx;
                a.y += dy;
            }
        }
        Some(m)
    }
}

/// Built-in ChiralCraft `Tools` ring stamps (`P15`, `P27`, `P29`, `P34b`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuiltinTemplate {
    Benzene,
    Cyclohexane,
    Cyclopentane,
    Cyclobutane,
    Cyclopropane,
    CyclohexaneChair,
    Bicyclo221Heptane,
}

fn mk_atom(id: u32, z: i32, x: f64, y: f64) -> ChemAtom {
    let mut a = ChemAtom::new(id, z);
    a.x = x;
    a.y = y;
    a
}

fn mk_bond(id: u32, begin_atom: usize, end_atom: usize, order: i32) -> ChemBond {
    ChemBond {
        id,
        begin_atom,
        end_atom,
        order,
        display: BondDisplay::Solid,
        visible: true,
        ignore_warnings: false,
    }
}

pub fn stamp_template(tpl: BuiltinTemplate, origin_x: f64, origin_y: f64, bond_len: f64) -> ChemMolecule {
    let mut m = ChemMolecule::default();
    let regular_poly = |n: usize, aromatic: bool| -> ChemMolecule {
        let mut p = ChemMolecule::default();
        let r = bond_len / (2.0 * (std::f64::consts::PI / (n as f64)).sin());
        for i in 0..n {
            let ang = (i as f64) * 2.0 * std::f64::consts::PI / (n as f64) - std::f64::consts::FRAC_PI_2;
            p.atoms
                .push(mk_atom((i + 1) as u32, 6, origin_x + r * ang.cos(), origin_y + r * ang.sin()));
        }
        for i in 0..n {
            let ord = if aromatic && i % 2 == 0 { 2 } else { 1 };
            p.bonds.push(mk_bond((i + 1) as u32, i, (i + 1) % n, ord));
        }
        p
    };
    match tpl {
        BuiltinTemplate::Benzene => return regular_poly(6, true),
        BuiltinTemplate::Cyclohexane => return regular_poly(6, false),
        BuiltinTemplate::Cyclopentane => return regular_poly(5, false),
        BuiltinTemplate::Cyclobutane => return regular_poly(4, false),
        BuiltinTemplate::Cyclopropane => return regular_poly(3, false),
        BuiltinTemplate::CyclohexaneChair => {
            let pts = [
                (0.0, 0.25),
                (0.95, 0.0),
                (1.90, 0.25),
                (1.65, 0.85),
                (0.70, 1.10),
                (-0.25, 0.85),
            ];
            for (i, &(px, py)) in pts.iter().enumerate() {
                m.atoms.push(mk_atom((i + 1) as u32, 6, origin_x + px * bond_len, origin_y + py * bond_len));
            }
            for i in 0..6 {
                m.bonds.push(mk_bond((i + 1) as u32, i, (i + 1) % 6, 1));
            }
        }
        BuiltinTemplate::Bicyclo221Heptane => {
            let pts = [
                (0.0, 0.5),
                (0.8, 0.0),
                (1.6, 0.5),
                (1.6, 1.4),
                (0.8, 1.8),
                (0.0, 1.4),
                (0.8, 0.9),
            ];
            for (i, &(px, py)) in pts.iter().enumerate() {
                m.atoms.push(mk_atom((i + 1) as u32, 6, origin_x + px * bond_len, origin_y + py * bond_len));
            }
            let edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (6, 3)];
            for (i, &(b, e)) in edges.iter().enumerate() {
                m.bonds.push(mk_bond((i + 1) as u32, b, e, 1));
            }
        }
    }
    m
}

// ============================================================================
// CDXML Extractors for Extended Scene Objects
// ============================================================================

fn parse_bbox_4(attr: Option<&str>) -> (f64, f64, f64, f64) {
    let Some(s) = attr else {
        return (0.0, 0.0, 0.0, 0.0);
    };
    let v: Vec<f64> = s.split_whitespace().filter_map(|t| t.parse().ok()).collect();
    if v.len() >= 4 {
        (v[0], v[1], v[2], v[3])
    } else {
        (0.0, 0.0, 0.0, 0.0)
    }
}

pub fn extract_arrows_from_cdxml(cdxml: &CdxmlDocument) -> Vec<SceneArrow> {
    let mut out = Vec::new();
    for node in cdxml.elements.iter().filter(|e| e.tag == "arrow") {
        let id = node.attr("id").and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        let parse_pt = |attr: &str| -> (f64, f64) {
            node.attr(attr)
                .map(|s| {
                    let v: Vec<f64> = s.split_whitespace().filter_map(|t| t.parse::<f64>().ok()).collect();
                    if v.len() >= 2 { (v[0], v[1]) } else { (0.0, 0.0) }
                })
                .unwrap_or((0.0, 0.0))
        };
        out.push(SceneArrow {
            id,
            tail: parse_pt("Tail3D"),
            head: parse_pt("Head3D"),
            arrowhead_head: node.attr("ArrowheadHead").unwrap_or("Full").to_string(),
            arrowhead_tail: node.attr("ArrowheadTail").map(|s| s.to_string()),
            no_go: node.attr("NoGo").map(|s| s.to_string()),
            dipole: node.attr("Dipole").map(|s| s == "yes").unwrap_or(false),
        });
    }
    out
}

pub fn extract_curves_from_cdxml(cdxml: &CdxmlDocument) -> Vec<SceneCurve> {
    let mut out = Vec::new();
    for node in cdxml.elements.iter().filter(|e| e.tag == "curve") {
        let id = node.attr("id").and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        let nums: Vec<f64> = node
            .attr("CurvePoints")
            .unwrap_or("")
            .split_whitespace()
            .filter_map(|t| t.parse().ok())
            .collect();
        let points: Vec<(f64, f64)> = nums.chunks_exact(2).map(|c| (c[0], c[1])).collect();
        let closed = node.attr("Closed").map(|s| s == "yes").unwrap_or(false);
        if !points.is_empty() {
            out.push(SceneCurve { id, points, closed });
        }
    }
    out
}

pub fn extract_graphics_from_cdxml(cdxml: &CdxmlDocument) -> Vec<SceneGraphic> {
    let mut out = Vec::new();
    for node in cdxml.elements.iter().filter(|e| e.tag == "graphic") {
        let gtype = node.attr("GraphicType").unwrap_or("Line");
        if gtype == "Oval" && node.attr("OvalType") == Some("Circle") {
            continue; // handled by `aromatic_circles`
        }
        let id = node.attr("id").and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        let sub_type = node
            .attr("OrbitalType")
            .or_else(|| node.attr("RectangleType"))
            .or_else(|| node.attr("SymbolType"))
            .or_else(|| node.attr("BracketType"))
            .map(|s| s.to_string());
        out.push(SceneGraphic {
            id,
            graphic_type: gtype.to_string(),
            sub_type,
            bbox: parse_bbox_4(node.attr("BoundingBox")),
        });
    }
    out
}

pub fn extract_tables_from_cdxml(cdxml: &CdxmlDocument) -> Vec<SceneTable> {
    let mut out = Vec::new();
    for node in cdxml.elements.iter().filter(|e| e.tag == "table") {
        let id = node.attr("id").and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        let pages: usize = node
            .children
            .iter()
            .filter(|c| matches!(c, CdxmlNodeChild::Element(cid) if cdxml.elements.get(*cid).map(|e| e.tag == "page").unwrap_or(false)))
            .count();
        let side = (pages as f64).sqrt().ceil().max(1.0) as usize;
        out.push(SceneTable {
            id,
            rows: side,
            cols: side,
            bbox: parse_bbox_4(node.attr("BoundingBox")),
            cells: vec![vec![String::new(); side]; side],
        });
    }
    out
}

pub fn extract_tlc_plates_from_cdxml(cdxml: &CdxmlDocument) -> Vec<SceneTlcPlate> {
    let mut out = Vec::new();
    for node in cdxml.elements.iter().filter(|e| e.tag == "tlcplate") {
        let id = node.attr("id").and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
        let origin_fraction = node
            .attr("OriginFraction")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.1);
        let mut spot_rfs = Vec::new();
        for c in &node.children {
            if let CdxmlNodeChild::Element(lane_id) = c {
                if let Some(lane) = cdxml.elements.get(*lane_id) {
                    if lane.tag == "tlclane" {
                        for sc in &lane.children {
                            if let CdxmlNodeChild::Element(sp_id) = sc {
                                if let Some(sp) = cdxml.elements.get(*sp_id) {
                                    if sp.tag == "tlcspot" {
                                        if let Some(rf) = sp.attr("Rf").and_then(|v| v.parse().ok()) {
                                            spot_rfs.push(rf);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let n_sp = spot_rfs.len();
        out.push(SceneTlcPlate {
            id,
            bbox: parse_bbox_4(node.attr("BoundingBox")),
            origin_fraction,
            spot_rfs,
            spot_labels: (1..=n_sp).map(|i| format!("L{i}")).collect(),
            spot_colors: vec![(25, 25, 25); n_sp],
            solvent_label: "Rf".to_string(),
        });
    }
    out
}

pub fn extract_page_header_footer(cdxml: &CdxmlDocument) -> PageHeaderFooter {
    if let Some(pg) = cdxml.elements.iter().find(|e| e.tag == "page") {
        PageHeaderFooter {
            header: pg.attr("Header").map(|s| s.to_string()),
            footer: pg.attr("Footer").map(|s| s.to_string()),
            width_pages: pg.attr("WidthPages").and_then(|s| s.parse().ok()).unwrap_or(1),
            height_pages: pg.attr("HeightPages").and_then(|s| s.parse().ok()).unwrap_or(1),
        }
    } else {
        PageHeaderFooter::default()
    }
}

// ============================================================================
// CDXML Serialization & Full Vector SVG Rendering for Extended Scene Objects
// ============================================================================

pub fn serialize_scene_ext_cdxml(
    arrows: &[SceneArrow],
    curves: &[SceneCurve],
    graphics: &[SceneGraphic],
    tlc_plates: &[SceneTlcPlate],
    out: &mut String,
) {
    let esc = |s: &str| -> String {
        s.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
    };
    for a in arrows {
        let tail_attr = a
            .arrowhead_tail
            .as_deref()
            .map(|t| format!(" ArrowheadTail=\"{}\"", esc(t)))
            .unwrap_or_default();
        let nogo_attr = a
            .no_go
            .as_deref()
            .map(|n| format!(" NoGo=\"{}\"", esc(n)))
            .unwrap_or_default();
        let dip_attr = if a.dipole { " Dipole=\"yes\"" } else { "" };
        out.push_str(&format!(
            "<arrow id=\"{}\" Tail3D=\"{:.2} {:.2} 0\" Head3D=\"{:.2} {:.2} 0\" ArrowheadType=\"Solid\" ArrowheadHead=\"{}\"{}{}{}/>",
            a.id, a.tail.0, a.tail.1, a.head.0, a.head.1, esc(&a.arrowhead_head), tail_attr, nogo_attr, dip_attr
        ));
    }
    for c in curves {
        let pts: Vec<String> = c.points.iter().map(|(x, y)| format!("{x:.2} {y:.2}")).collect();
        let closed = if c.closed { " Closed=\"yes\"" } else { "" };
        out.push_str(&format!(
            "<curve id=\"{}\" CurvePoints=\"{}\"{}/>",
            c.id,
            pts.join(" "),
            closed
        ));
    }
    for g in graphics {
        let sub = match (g.graphic_type.as_str(), g.sub_type.as_deref()) {
            ("Orbital", Some(s)) => format!(" OrbitalType=\"{}\"", esc(s)),
            ("Rectangle", Some(s)) => format!(" RectangleType=\"{}\"", esc(s)),
            ("Symbol", Some(s)) => format!(" SymbolType=\"{}\"", esc(s)),
            ("Bracket", Some(s)) => format!(" BracketType=\"{}\"", esc(s)),
            _ => String::new(),
        };
        out.push_str(&format!(
            "<graphic id=\"{}\" GraphicType=\"{}\"{} BoundingBox=\"{:.2} {:.2} {:.2} {:.2}\"/>",
            g.id, esc(&g.graphic_type), sub, g.bbox.0, g.bbox.1, g.bbox.2, g.bbox.3
        ));
    }
    for tlc in tlc_plates {
        out.push_str(&format!(
            "<tlcplate id=\"{}\" BoundingBox=\"{:.2} {:.2} {:.2} {:.2}\" OriginFraction=\"{:.2}\">",
            tlc.id, tlc.bbox.0, tlc.bbox.1, tlc.bbox.2, tlc.bbox.3, tlc.origin_fraction
        ));
        if !tlc.spot_rfs.is_empty() {
            out.push_str("<tlclane>");
            for rf in &tlc.spot_rfs {
                out.push_str(&format!("<tlcspot Rf=\"{rf:.2}\"/>"));
            }
            out.push_str("</tlclane>");
        }
        out.push_str("</tlcplate>");
    }
}

/// Full vector SVG renderer for arrows (with solid/half/tail arrowheads & `NoGo="Cross"`),
/// Bezier curves, graphics (orbitals/rectangles/symbols), tables, TLC plates, and page headers/footers.
pub fn arrows_to_svg_fragment(arrows: &[SceneArrow]) -> String {
    arrows_to_svg_fragment_styled(arrows, "black")
}

pub fn arrows_to_svg_fragment_styled(arrows: &[SceneArrow], ink: &str) -> String {
    let mut out = String::new();
    for a in arrows {
        let (x1, y1) = a.tail;
        let (x2, y2) = a.head;
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = dx.hypot(dy).max(1e-6);
        let ux = dx / len;
        let uy = dy / len;
        let nx = -uy;
        let ny = ux;
        let hl = 5.0;
        let hw = 2.0;

        // Equilibrium double-shaft if both head & tail are HalfLeft
        let is_eq = a.arrowhead_head == "HalfLeft" && a.arrowhead_tail.as_deref() == Some("HalfLeft");
        if is_eq {
            let off = 1.6;
            out.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n",
                x1 + nx * off, y1 + ny * off, x2 + nx * off, y2 + ny * off
            ));
            out.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n",
                x1 - nx * off, y1 - ny * off, x2 - nx * off, y2 - ny * off
            ));
        } else {
            out.push_str(&format!(
                "<line x1=\"{x1:.2}\" y1=\"{y1:.2}\" x2=\"{x2:.2}\" y2=\"{y2:.2}\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
            ));
        }

        // Head arrowhead
        if a.arrowhead_head == "HalfLeft" {
            out.push_str(&format!(
                "<polygon points=\"{x2:.2},{y2:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{ink}\"/>\n",
                x2 - ux * hl + nx * hw,
                y2 - uy * hl + ny * hw,
                x2 - ux * (hl * 0.7),
                y2 - uy * (hl * 0.7)
            ));
        } else {
            out.push_str(&format!(
                "<polygon points=\"{x2:.2},{y2:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{ink}\"/>\n",
                x2 - ux * hl + nx * hw,
                y2 - uy * hl + ny * hw,
                x2 - ux * hl - nx * hw,
                y2 - uy * hl - ny * hw
            ));
        }

        // Optional tail arrowhead (resonance / equilibrium)
        if let Some(tail_kind) = a.arrowhead_tail.as_deref() {
            if tail_kind == "HalfLeft" {
                out.push_str(&format!(
                    "<polygon points=\"{x1:.2},{y1:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{ink}\"/>\n",
                    x1 + ux * hl - nx * hw,
                    y1 + uy * hl - ny * hw,
                    x1 + ux * (hl * 0.7),
                    y1 + uy * (hl * 0.7)
                ));
            } else {
                out.push_str(&format!(
                    "<polygon points=\"{x1:.2},{y1:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{ink}\"/>\n",
                    x1 + ux * hl + nx * hw,
                    y1 + uy * hl + ny * hw,
                    x1 + ux * hl - nx * hw,
                    y1 + uy * hl - ny * hw
                ));
            }
        }

        // NoGo="Cross" X mark across midpoint
        if a.no_go.as_deref() == Some("Cross") {
            let mx = 0.5 * (x1 + x2);
            let my = 0.5 * (y1 + y2);
            out.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"0.75\"/>\n",
                mx - 4.0, my - 4.0, mx + 4.0, my + 4.0
            ));
            out.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"0.75\"/>\n",
                mx - 4.0, my + 4.0, mx + 4.0, my - 4.0
            ));
        }
    }
    out
}

pub fn render_scene_ext_svg(
    arrows: &[SceneArrow],
    curves: &[SceneCurve],
    graphics: &[SceneGraphic],
    tables: &[SceneTable],
    tlc_plates: &[SceneTlcPlate],
    hf: &PageHeaderFooter,
) -> String {
    render_scene_ext_svg_styled(arrows, curves, graphics, tables, tlc_plates, hf, "black")
}

pub fn render_scene_ext_svg_styled(
    arrows: &[SceneArrow],
    curves: &[SceneCurve],
    graphics: &[SceneGraphic],
    tables: &[SceneTable],
    tlc_plates: &[SceneTlcPlate],
    hf: &PageHeaderFooter,
    ink: &str,
) -> String {
    let mut svg = arrows_to_svg_fragment_styled(arrows, ink);
    for c in curves {
        if c.points.len() >= 2 {
            let mut d = format!("M {:.2} {:.2}", c.points[0].0, c.points[0].1);
            let mut i = 1;
            while i + 2 < c.points.len() {
                d.push_str(&format!(
                    " C {:.2} {:.2}, {:.2} {:.2}, {:.2} {:.2}",
                    c.points[i].0,
                    c.points[i].1,
                    c.points[i + 1].0,
                    c.points[i + 1].1,
                    c.points[i + 2].0,
                    c.points[i + 2].1
                ));
                i += 3;
            }
            while i < c.points.len() {
                d.push_str(&format!(" L {:.2} {:.2}", c.points[i].0, c.points[i].1));
                i += 1;
            }
            if c.closed {
                d.push_str(" Z");
            }
            svg.push_str(&format!(
                "<path d=\"{d}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
            ));
        }
    }
    for g in graphics {
        let (x1, y1, x2, y2) = g.bbox;
        let w = (x2 - x1).abs();
        let h = (y2 - y1).abs();
        let rx = x1.min(x2);
        let ry = y1.min(y2);
        if g.graphic_type == "Rectangle" {
            svg.push_str(&format!(
                "<rect x=\"{rx:.2}\" y=\"{ry:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
            ));
        } else if g.graphic_type == "Orbital" || g.graphic_type == "Oval" {
            svg.push_str(&format!(
                "<ellipse cx=\"{:.2}\" cy=\"{:.2}\" rx=\"{:.2}\" ry=\"{:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n",
                rx + w * 0.5, ry + h * 0.5, (w * 0.5).max(1.0), (h * 0.5).max(1.0)
            ));
        } else if g.graphic_type == "Symbol" {
            let cx = rx + w * 0.5;
            let cy = ry + h * 0.5;
            let r = (w.min(h) * 0.5).max(4.0);
            svg.push_str(&format!(
                "<circle cx=\"{cx:.2}\" cy=\"{cy:.2}\" r=\"{r:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
            ));
        }
    }
    for t in tables {
        let (x1, y1, x2, y2) = t.bbox;
        let rx = x1.min(x2);
        let ry = y1.min(y2);
        let w = (x2 - x1).abs().max(24.0);
        let h = (y2 - y1).abs().max(18.0);
        let nr = t.rows.max(1);
        let nc = t.cols.max(1);
        svg.push_str(&format!(
            "<rect x=\"{rx:.2}\" y=\"{ry:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
        ));
        for r in 1..nr {
            let y = ry + h * (r as f64) / (nr as f64);
            svg.push_str(&format!("<line x1=\"{rx:.2}\" y1=\"{y:.2}\" x2=\"{:.2}\" y2=\"{y:.2}\" stroke=\"{ink}\" stroke-width=\"0.50\"/>\n", rx + w));
        }
        for c in 1..nc {
            let x = rx + w * (c as f64) / (nc as f64);
            svg.push_str(&format!("<line x1=\"{x:.2}\" y1=\"{ry:.2}\" x2=\"{x:.2}\" y2=\"{:.2}\" stroke=\"{ink}\" stroke-width=\"0.50\"/>\n", ry + h));
        }
        for (ri, row) in t.cells.iter().take(nr).enumerate() {
            for (ci, txt) in row.iter().take(nc).enumerate() {
                if !txt.is_empty() {
                    let tx = rx + w * ((ci as f64) + 0.5) / (nc as f64);
                    let ty = ry + h * ((ri as f64) + 0.5) / (nr as f64);
                    let esc = txt.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
                    svg.push_str(&format!(
                        "<text x=\"{tx:.2}\" y=\"{ty:.2}\" font-size=\"8\" text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"{ink}\">{esc}</text>\n"
                    ));
                }
            }
        }
    }
    let muted_ink = if ink == "black" { "#444444" } else { ink };
    for tlc in tlc_plates {
        let (x1, y1, x2, y2) = tlc.bbox;
        let rx = x1.min(x2);
        let ry = y1.min(y2);
        let w = (x2 - x1).abs().max(20.0);
        let h = (y2 - y1).abs().max(40.0);
        let oy = ry + h * (1.0 - tlc.origin_fraction);
        svg.push_str(&format!(
            "<rect x=\"{rx:.2}\" y=\"{ry:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"none\" stroke=\"{ink}\" stroke-width=\"0.60\"/>\n"
        ));
        svg.push_str(&format!(
            "<line x1=\"{rx:.2}\" y1=\"{oy:.2}\" x2=\"{:.2}\" y2=\"{oy:.2}\" stroke=\"{ink}\" stroke-dasharray=\"2,2\" stroke-width=\"0.50\"/>\n",
            rx + w
        ));
        for (i, &rf) in tlc.spot_rfs.iter().enumerate() {
            let sx = rx + w * ((i + 1) as f64) / ((tlc.spot_rfs.len() + 1) as f64);
            let sy = oy - rf * (oy - ry);
            let (cr, cg, cb) = tlc.spot_colors.get(i).copied().unwrap_or((25, 25, 25));
            svg.push_str(&format!(
                "<circle cx=\"{sx:.2}\" cy=\"{sy:.2}\" r=\"2.4\" fill=\"rgb({cr},{cg},{cb})\"/>\n"
            ));
            if let Some(lbl) = tlc.spot_labels.get(i) {
                if !lbl.is_empty() {
                    let esc = lbl.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
                    svg.push_str(&format!(
                        "<text x=\"{sx:.2}\" y=\"{:.2}\" font-size=\"6.5\" text-anchor=\"middle\" fill=\"{ink}\">{esc}</text>\n",
                        oy + 7.5
                    ));
                }
            }
        }
        if !tlc.solvent_label.is_empty() {
            let esc = tlc.solvent_label.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
            svg.push_str(&format!(
                "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"6.5\" text-anchor=\"middle\" fill=\"{muted_ink}\">{esc}</text>\n",
                rx + w * 0.5,
                ry - 3.0
            ));
        }
    }
    if let Some(hdr) = &hf.header {
        let esc = hdr.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
        svg.push_str(&format!(
            "<text x=\"270.00\" y=\"24.00\" font-size=\"9\" text-anchor=\"middle\" fill=\"{ink}\">{esc}</text>\n"
        ));
    }
    if let Some(ftr) = &hf.footer {
        let esc = ftr.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
        svg.push_str(&format!(
            "<text x=\"270.00\" y=\"700.00\" font-size=\"9\" text-anchor=\"middle\" fill=\"{ink}\">{esc}</text>\n"
        ));
    }
    svg
}

// ============================================================================
// Pure-Geometry Special Chemical Cores (Haworth Sugars, Macrocycles, Bridged Cages)
// ============================================================================

/// Builds a textbook Haworth projection of β-D-Glucopyranose with bold/wedged front edge & vertical substituents.
pub fn build_haworth_pyranose(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let s = blen;
    // Perspective flat hexagon: C1(right), C2(front-right), C3(front-left), C4(left), C5(back-left), O5(back-right)
    let coords = [
        (1.00 * s, 0.00 * s, 6),  // 0: C1
        (0.55 * s, 0.55 * s, 6),  // 1: C2
        (-0.55 * s, 0.55 * s, 6), // 2: C3
        (-1.00 * s, 0.00 * s, 6), // 3: C4
        (-0.55 * s, -0.48 * s, 6),// 4: C5
        (0.55 * s, -0.48 * s, 8), // 5: O5
        // Vertical axial substituents (up/down)
        (1.00 * s, -0.82 * s, 8), // 6: C1-OH (β, up)
        (0.55 * s, 1.38 * s, 8),  // 7: C2-OH (down)
        (-0.55 * s, -0.28 * s, 8),// 8: C3-OH (up)
        (-1.00 * s, 0.82 * s, 8), // 9: C4-OH (down)
        (-0.55 * s, -1.25 * s, 6),// 10: C6 (up)
        (-1.25 * s, -1.60 * s, 8),// 11: C6-OH
    ];
    for (i, &(x, y, z)) in coords.iter().enumerate() {
        let mut a = ChemAtom::new((i + 1) as u32, z);
        a.x = x;
        a.y = y;
        m.add_atom(a);
    }
    let bonds = [
        (0, 1, 1, BondDisplay::WedgeBegin),
        (1, 2, 1, BondDisplay::Bold),
        (2, 3, 1, BondDisplay::WedgeEnd),
        (3, 4, 1, BondDisplay::Solid),
        (4, 5, 1, BondDisplay::Solid),
        (5, 0, 1, BondDisplay::Solid),
        (0, 6, 1, BondDisplay::Solid),
        (1, 7, 1, BondDisplay::Solid),
        (2, 8, 1, BondDisplay::Solid),
        (3, 9, 1, BondDisplay::Solid),
        (4, 10, 1, BondDisplay::Solid),
        (10, 11, 1, BondDisplay::Solid),
    ];
    for (k, &(a, b, ord, disp)) in bonds.iter().enumerate() {
        m.add_bond((20 + k) as u32, a, b, ord, disp);
    }
    m
}

/// Builds a textbook Haworth projection of β-D-Ribofuranose with bold front edge & vertical substituents.
pub fn build_haworth_furanose(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let s = blen;
    let coords = [
        (0.00 * s, -0.62 * s, 8),  // 0: O4 (top)
        (0.95 * s, -0.05 * s, 6),  // 1: C1 (right)
        (0.55 * s, 0.52 * s, 6),   // 2: C2 (front-right)
        (-0.55 * s, 0.52 * s, 6),  // 3: C3 (front-left)
        (-0.95 * s, -0.05 * s, 6), // 4: C4 (left)
        (0.95 * s, -0.88 * s, 8),  // 5: C1-OH (up)
        (0.55 * s, 1.32 * s, 8),   // 6: C2-OH (down)
        (-0.55 * s, 1.32 * s, 8),  // 7: C3-OH (down)
        (-0.95 * s, -0.88 * s, 6), // 8: C5 (up)
        (-1.62 * s, -1.25 * s, 8), // 9: C5-OH
    ];
    for (i, &(x, y, z)) in coords.iter().enumerate() {
        let mut a = ChemAtom::new((i + 1) as u32, z);
        a.x = x;
        a.y = y;
        m.add_atom(a);
    }
    let bonds = [
        (0, 1, 1, BondDisplay::Solid),
        (1, 2, 1, BondDisplay::WedgeBegin),
        (2, 3, 1, BondDisplay::Bold),
        (3, 4, 1, BondDisplay::WedgeEnd),
        (4, 0, 1, BondDisplay::Solid),
        (1, 5, 1, BondDisplay::Solid),
        (2, 6, 1, BondDisplay::Solid),
        (3, 7, 1, BondDisplay::Solid),
        (4, 8, 1, BondDisplay::Solid),
        (8, 9, 1, BondDisplay::Solid),
    ];
    for (k, &(a, b, ord, disp)) in bonds.iter().enumerate() {
        m.add_bond((20 + k) as u32, a, b, ord, disp);
    }
    m
}

/// Builds a symmetric planar 24-atom Porphine (Porphyrin macrocycle) core!
pub fn build_porphine_core(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let s = blen;
    // Four pyrrole quadrants rotated by 0°, 90°, 180°, 270° + four meso CH bridges at 45°, 135°, 225°, 315°
    let mut quad_a = [0usize; 4];
    let mut quad_b = [0usize; 4];
    let mut next_aid = 1u32;
    let mut next_bid = 100u32;

    for q in 0..4 {
        let ang = (q as f64) * std::f64::consts::FRAC_PI_2;
        let (sa, ca) = ang.sin_cos();
        let rot = |lx: f64, ly: f64| -> (f64, f64) { (lx * ca - ly * sa, lx * sa + ly * ca) };
        // Local pyrrole facing inward toward origin: N at (1.15s, 0), alpha C at (1.75s, ±0.78s), beta C at (2.75s, ±0.48s)
        let pts = [
            (1.12 * s, 0.0, 7),
            (1.78 * s, -0.78 * s, 6),
            (2.75 * s, -0.48 * s, 6),
            (2.75 * s, 0.48 * s, 6),
            (1.78 * s, 0.78 * s, 6),
        ];
        let base = m.atoms.len();
        for &(lx, ly, z) in &pts {
            let (rx, ry) = rot(lx, ly);
            let mut a = ChemAtom::new(next_aid, z);
            next_aid += 1;
            a.x = rx;
            a.y = ry;
            m.add_atom(a);
        }
        quad_a[q] = base + 1;
        quad_b[q] = base + 4;
        // Valid 11-double-bond 18π Kekulé assignment for free-base porphine (C20H14N4, 0 valence warnings):
        // q=0 (NH): 1=2, 3=4 | q=1 (N): 2=3, 4=0 | q=2 (NH): 2=3 | q=3 (N): 0=1, 2=3
        let ring_b = match q {
            0 => [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 0, 1)],
            1 => [(0, 1, 1), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 0, 2)],
            2 => [(0, 1, 1), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 0, 1)],
            _ => [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 0, 1)],
        };
        for &(u, v, ord) in &ring_b {
            m.add_bond(next_bid, base + u, base + v, ord, BondDisplay::Solid);
            next_bid += 1;
        }
    }
    // Meso methine bridges: meso0=a1, meso1=a2, b2=meso2, b3=meso3
    for q in 0..4 {
        let ang = (q as f64) * std::f64::consts::FRAC_PI_2 + std::f64::consts::FRAC_PI_4;
        let mut meso = ChemAtom::new(next_aid, 6);
        next_aid += 1;
        meso.x = 2.02 * s * ang.cos();
        meso.y = 2.02 * s * ang.sin();
        let mi = m.add_atom(meso);
        let (ord_in, ord_out) = if q < 2 { (1, 2) } else { (2, 1) };
        m.add_bond(next_bid, quad_b[q], mi, ord_in, BondDisplay::Solid);
        next_bid += 1;
        m.add_bond(next_bid, mi, quad_a[(q + 1) % 4], ord_out, BondDisplay::Solid);
        next_bid += 1;
    }
    m
}

/// Builds a symmetric 18-Crown-6 macrocycle (`(-CH2-CH2-O-)6`).
pub fn build_crown_ether_18_6(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let n = 18usize;
    let r = blen / (2.0 * (std::f64::consts::PI / (n as f64)).sin());
    for i in 0..n {
        let ang = (i as f64) * std::f64::consts::TAU / (n as f64) - std::f64::consts::FRAC_PI_2;
        // Slight inward puckering on the 6 oxygens for textbook crown look
        let ri = if i % 3 == 0 { r * 0.88 } else { r };
        let z = if i % 3 == 0 { 8 } else { 6 };
        let mut a = ChemAtom::new((i + 1) as u32, z);
        a.x = ri * ang.cos();
        a.y = ri * ang.sin();
        m.add_atom(a);
    }
    for i in 0..n {
        m.add_bond((30 + i) as u32, i, (i + 1) % n, 1, BondDisplay::Solid);
    }
    m
}

/// Builds a classic 3D-projected Adamantane (`C10H16`) diamondoid cage with true tetrahedral
/// symmetry, Z-depth coordinates, and textbook crossing-bond bridge clearance!
pub fn build_adamantane_2d(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let scale = blen / 3.0_f64.sqrt();
    // Exact tetrahedral diamond-lattice coordinates for C10H16 (all 12 bonds = sqrt(3))
    let raw: [(f64, f64, f64); 10] = [
        (1.0, 1.0, 1.0),    // 0: bridgehead
        (1.0, -1.0, -1.0),  // 1: bridgehead
        (-1.0, 1.0, -1.0),  // 2: bridgehead
        (-1.0, -1.0, 1.0),  // 3: bridgehead
        (2.0, 0.0, 0.0),    // 4: +x CH2 (0-1)
        (-2.0, 0.0, 0.0),   // 5: -x CH2 (2-3)
        (0.0, 2.0, 0.0),    // 6: +y CH2 (0-2)
        (0.0, -2.0, 0.0),   // 7: -y CH2 (1-3)
        (0.0, 0.0, 2.0),    // 8: +z CH2 (0-3)
        (0.0, 0.0, -2.0),   // 9: -z CH2 (1-2)
    ];
    let yaw = 0.48_f64;   // ~27.5 deg
    let pitch = 0.36_f64; // ~20.6 deg
    let (sy, cy) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    for (i, &(rx, ry, rz)) in raw.iter().enumerate() {
        let x1 = rx * cy + rz * sy;
        let z1 = -rx * sy + rz * cy;
        let y2 = ry * cp - z1 * sp;
        let z2 = ry * sp + z1 * cp;
        let mut a = ChemAtom::new((i + 1) as u32, 6);
        a.x = x1 * scale;
        a.y = -y2 * scale;
        a.z_coord = z2 * scale;
        m.add_atom(a);
    }
    m.has_3d_coords = true;
    let edges = [
        (0, 4), (4, 1),
        (2, 5), (5, 3),
        (0, 6), (6, 2),
        (1, 7), (7, 3),
        (0, 8), (8, 3),
        (1, 9), (9, 2),
    ];
    for (k, &(u, v)) in edges.iter().enumerate() {
        m.add_bond((20 + k) as u32, u, v, 1, BondDisplay::Solid);
    }
    m
}

/// Builds a classic 2D bridged Norbornane (`bicyclo[2.2.1]heptane`) with clean crossing-bond geometry!
pub fn build_norbornane_2d(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let s = blen;
    let pts = [
        (-0.90 * s, 0.20 * s, 0.2),   // 0: left bridgehead
        (-0.52 * s, 0.85 * s, 0.6),   // 1: front-left CH2
        (0.52 * s, 0.85 * s, 0.6),    // 2: front-right CH2
        (0.90 * s, 0.20 * s, 0.2),    // 3: right bridgehead
        (0.60 * s, -0.30 * s, -0.5),  // 4: back-right CH2
        (-0.60 * s, -0.30 * s, -0.5), // 5: back-left CH2
        (0.00 * s, -0.82 * s, 0.9),   // 6: C7 apex bridge
    ];
    for (i, &(x, y, zc)) in pts.iter().enumerate() {
        let mut a = ChemAtom::new((i + 1) as u32, 6);
        a.x = x;
        a.y = y;
        a.z_coord = zc * s;
        m.add_atom(a);
    }
    m.has_3d_coords = true;
    let edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (6, 3)];
    for (k, &(u, v)) in edges.iter().enumerate() {
        m.add_bond((20 + k) as u32, u, v, 1, BondDisplay::Solid);
    }
    m
}
