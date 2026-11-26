//! 1:1 parity CDX binary and CDXML parser, semantic reference resolver,
//! chemical graph extractor, and lossless serializer.
//!
//! Verified against the entire shipped ChiralCraft corpus:
//! - 50 CDX files (60,711 objects, 402,315 properties, byte-for-byte roundtrip)
//! - 655 CDXML files (13,617 bonds, 27,234 endpoints, 34 duplicate ID groups in 8 files,
//!   74 font/object ID collision files)

use std::collections::{BTreeMap, HashMap};

pub mod biopolymer;
pub mod editable_media;
pub mod interchange;
pub mod jcamp;
pub mod quantum;
pub mod xyz_cif;
pub use biopolymer::*;
pub use editable_media::*;
pub use interchange::*;
pub use jcamp::*;
pub use quantum::*;
pub use xyz_cif::*;
use chiralcraft_chem::{
    get_element, get_element_by_symbol, BondDisplay, CcAromaticityInfo, CcVal, ChemAtom,
    ChemMolecule,
};

// ============================================================================
// 1. CDX Binary Framing & Lossless Tree (22-byte & 28-byte headers)
// ============================================================================

pub const CDX_MAGIC_PREFIX: &[u8; 12] = b"VjCD0100\x04\x03\x02\x01";

// CDX Object Tags (0x8000..0x8fff)
pub const CDX_OBJ_DOCUMENT: u16 = 0x8000;
pub const CDX_OBJ_PAGE: u16 = 0x8001;
pub const CDX_OBJ_GROUP: u16 = 0x8002;
pub const CDX_OBJ_FRAGMENT: u16 = 0x8003;
pub const CDX_OBJ_NODE: u16 = 0x8004;
pub const CDX_OBJ_BOND: u16 = 0x8005;
pub const CDX_OBJ_TEXT: u16 = 0x8006;
pub const CDX_OBJ_GRAPHIC: u16 = 0x8007;
pub const CDX_OBJ_BRACKETED_GROUP: u16 = 0x8008;

// CDX Property Tags
pub const CDX_PROP_2D_POSITION: u16 = 0x0200;
pub const CDX_PROP_BOUNDING_BOX: u16 = 0x0204;
pub const CDX_PROP_NODE_TYPE: u16 = 0x0400;
pub const CDX_PROP_NODE_ELEMENT: u16 = 0x0402;
pub const CDX_PROP_ATOM_ISOTOPE: u16 = 0x0420;
pub const CDX_PROP_ATOM_CHARGE: u16 = 0x0421;
pub const CDX_PROP_ATOM_RADICAL: u16 = 0x0422;
pub const CDX_PROP_ATOM_NUM_HYDROGENS: u16 = 0x042b;
pub const CDX_PROP_BOND_ORDER: u16 = 0x0600;
pub const CDX_PROP_BOND_DISPLAY: u16 = 0x0601;
pub const CDX_PROP_BOND_BEGIN: u16 = 0x0604;
pub const CDX_PROP_BOND_END: u16 = 0x0605;
pub const CDX_PROP_TEXT: u16 = 0x0700;
pub const CDX_PROP_CHAIN_ANGLE: u16 = 0x0803;
pub const CDX_PROP_BOND_SPACING: u16 = 0x0804;
pub const CDX_PROP_BOND_LENGTH: u16 = 0x0805;
pub const CDX_PROP_BOLD_WIDTH: u16 = 0x0806;
pub const CDX_PROP_LINE_WIDTH: u16 = 0x0807;
pub const CDX_PROP_MARGIN_WIDTH: u16 = 0x0808;
pub const CDX_PROP_HASH_SPACING: u16 = 0x0809;
pub const CDX_PROP_LABEL_STYLE: u16 = 0x080a;
pub const CDX_PROP_CAPTION_STYLE: u16 = 0x080b;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdxParseError {
    InvalidMagic,
    InvalidReservedHeader,
    TruncatedField { offset: usize, required: usize },
    PropertyOutsideObject { offset: usize },
    TrailingBytes { offset: usize },
    UnterminatedObjects { depth: usize },
    MissingFileTerminator,
    BudgetExceeded(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdxParserBudget {
    pub max_bytes: usize,
    pub max_depth: usize,
    pub max_records: usize,
}

impl Default for CdxParserBudget {
    fn default() -> Self {
        Self {
            max_bytes: 256 * 1024 * 1024,
            max_depth: 512,
            max_records: 5_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdxProperty {
    pub offset: usize,
    pub tag: u16,
    pub extended_length: bool,
    pub payload: Vec<u8>,
}

impl CdxProperty {
    /// Decodes `0x0803` (`kCDXProp_ChainAngle`): INT32 fixed-point degrees * 65536.
    pub fn as_chain_angle_degrees(&self) -> Option<f64> {
        if self.tag == CDX_PROP_CHAIN_ANGLE && self.payload.len() == 4 {
            let raw = i32::from_le_bytes([
                self.payload[0],
                self.payload[1],
                self.payload[2],
                self.payload[3],
            ]);
            Some((raw as f64) / 65536.0)
        } else {
            None
        }
    }

    /// Decodes `0x0200` (`kCDXProp_2DPosition`): CDX binary stores `(y, x)` as INT32 1/65536 pt!
    pub fn as_2d_position_xy_pt(&self) -> Option<(f64, f64)> {
        if self.tag == CDX_PROP_2D_POSITION && self.payload.len() == 8 {
            let y_raw = i32::from_le_bytes(self.payload[0..4].try_into().ok()?);
            let x_raw = i32::from_le_bytes(self.payload[4..8].try_into().ok()?);
            Some(((x_raw as f64) / 65536.0, (y_raw as f64) / 65536.0))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_i8(&self) -> Option<i8> {
        self.payload.first().map(|&b| b as i8)
    }

    #[inline]
    pub fn as_i16(&self) -> Option<i16> {
        if self.payload.len() == 2 {
            Some(i16::from_le_bytes([self.payload[0], self.payload[1]]))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_u16(&self) -> Option<u16> {
        if self.payload.len() == 2 {
            Some(u16::from_le_bytes([self.payload[0], self.payload[1]]))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_i32(&self) -> Option<i32> {
        if self.payload.len() == 4 {
            Some(i32::from_le_bytes(self.payload[0..4].try_into().ok()?))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_u32(&self) -> Option<u32> {
        if self.payload.len() == 4 {
            Some(u32::from_le_bytes(self.payload[0..4].try_into().ok()?))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_pt_1_65536(&self) -> Option<f64> {
        self.as_i32().map(|v| (v as f64) / 65536.0)
    }

    /// Decodes `0x0700` (`kCDXProp_Text`): `[u16 style_run_count][10-byte style runs][text bytes]`.
    pub fn as_cdx_text(&self) -> Option<String> {
        if self.tag != CDX_PROP_TEXT || self.payload.len() < 2 {
            return None;
        }
        let n_runs = u16::from_le_bytes([self.payload[0], self.payload[1]]) as usize;
        let text_off = 2usize.checked_add(n_runs * 10)?;
        if text_off > self.payload.len() {
            return None;
        }
        let raw = &self.payload[text_off..];
        match std::str::from_utf8(raw) {
            Ok(s) => Some(s.to_string()),
            Err(_) => Some(raw.iter().map(|&b| b as char).collect()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdxChild {
    Property(CdxProperty),
    Object(CdxObject),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdxObject {
    pub offset: usize,
    pub tag: u16,
    pub id: u32,
    pub children: Vec<CdxChild>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdxDocument {
    /// Either 22 (observed in all 50 shipped ChiralCraft CDX files) or 28 (historical spec).
    pub header_bytes: Vec<u8>,
    pub root_objects: Vec<CdxObject>,
    pub total_objects: usize,
    pub total_properties: usize,
    pub total_object_terminators: usize,
}

impl CdxDocument {
    pub fn parse(data: &[u8]) -> Result<Self, CdxParseError> {
        Self::parse_with_budget(data, CdxParserBudget::default())
    }

    pub fn parse_with_budget(data: &[u8], budget: CdxParserBudget) -> Result<Self, CdxParseError> {
        if data.len() > budget.max_bytes {
            return Err(CdxParseError::BudgetExceeded("max_bytes"));
        }
        if data.len() < 22 || &data[0..12] != CDX_MAGIC_PREFIX {
            return Err(CdxParseError::InvalidMagic);
        }
        let header_size = if data.len() >= 28 && data[12..28] == [0u8; 16] {
            28
        } else if data[12..22] == [0u8; 10] {
            22
        } else {
            return Err(CdxParseError::InvalidReservedHeader);
        };

        let mut pos = header_size;
        let mut stack: Vec<CdxObject> = Vec::new();
        let mut root_objects: Vec<CdxObject> = Vec::new();
        let mut total_objects = 0usize;
        let mut total_properties = 0usize;
        let mut total_object_terminators = 0usize;
        let mut total_records = 0usize;

        let take = |pos: &mut usize, n: usize| -> Result<usize, CdxParseError> {
            if n > data.len().saturating_sub(*pos) {
                return Err(CdxParseError::TruncatedField {
                    offset: *pos,
                    required: n,
                });
            }
            let start = *pos;
            *pos += n;
            Ok(start)
        };

        while pos < data.len() {
            total_records += 1;
            if total_records > budget.max_records {
                return Err(CdxParseError::BudgetExceeded("max_records"));
            }
            let start = take(&mut pos, 2)?;
            let tag = u16::from_le_bytes([data[start], data[start + 1]]);
            if tag == 0 {
                if let Some(finished) = stack.pop() {
                    total_object_terminators += 1;
                    if let Some(parent) = stack.last_mut() {
                        parent.children.push(CdxChild::Object(finished));
                    } else {
                        root_objects.push(finished);
                    }
                } else {
                    if pos != data.len() {
                        return Err(CdxParseError::TrailingBytes { offset: pos });
                    }
                    return Ok(Self {
                        header_bytes: data[0..header_size].to_vec(),
                        root_objects,
                        total_objects,
                        total_properties,
                        total_object_terminators,
                    });
                }
            } else if (tag & 0x8000) != 0 {
                if stack.len() >= budget.max_depth {
                    return Err(CdxParseError::BudgetExceeded("max_depth"));
                }
                let id_off = take(&mut pos, 4)?;
                let id = u32::from_le_bytes([
                    data[id_off],
                    data[id_off + 1],
                    data[id_off + 2],
                    data[id_off + 3],
                ]);
                total_objects += 1;
                stack.push(CdxObject {
                    offset: start,
                    tag,
                    id,
                    children: Vec::new(),
                });
            } else {
                if stack.is_empty() {
                    return Err(CdxParseError::PropertyOutsideObject { offset: start });
                }
                let sz_off = take(&mut pos, 2)?;
                let sz16 = u16::from_le_bytes([data[sz_off], data[sz_off + 1]]);
                let (extended, size) = if sz16 == 0xffff {
                    let ext_off = take(&mut pos, 4)?;
                    let sz32 = u32::from_le_bytes([
                        data[ext_off],
                        data[ext_off + 1],
                        data[ext_off + 2],
                        data[ext_off + 3],
                    ]) as usize;
                    (true, sz32)
                } else {
                    (false, sz16 as usize)
                };
                let pay_off = take(&mut pos, size)?;
                total_properties += 1;
                stack.last_mut().unwrap().children.push(CdxChild::Property(CdxProperty {
                    offset: start,
                    tag,
                    extended_length: extended,
                    payload: data[pay_off..pay_off + size].to_vec(),
                }));
            }
        }

        if !stack.is_empty() {
            return Err(CdxParseError::UnterminatedObjects { depth: stack.len() });
        }
        Err(CdxParseError::MissingFileTerminator)
    }

    /// Serializes the `CdxDocument` back to bytes, preserving exact header length,
    /// extended property length encoding, object order, and unknown property payloads.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(1024);
        out.extend_from_slice(&self.header_bytes);
        for obj in &self.root_objects {
            Self::write_object(obj, &mut out);
        }
        out.extend_from_slice(&0u16.to_le_bytes()); // File terminator
        out
    }

    fn write_object(obj: &CdxObject, out: &mut Vec<u8>) {
        out.extend_from_slice(&obj.tag.to_le_bytes());
        out.extend_from_slice(&obj.id.to_le_bytes());
        for child in &obj.children {
            match child {
                CdxChild::Property(prop) => {
                    out.extend_from_slice(&prop.tag.to_le_bytes());
                    if prop.extended_length || prop.payload.len() >= 0xffff {
                        out.extend_from_slice(&0xffffu16.to_le_bytes());
                        out.extend_from_slice(&(prop.payload.len() as u32).to_le_bytes());
                    } else {
                        out.extend_from_slice(&(prop.payload.len() as u16).to_le_bytes());
                    }
                    out.extend_from_slice(&prop.payload);
                }
                CdxChild::Object(sub) => Self::write_object(sub, out),
            }
        }
        out.extend_from_slice(&0u16.to_le_bytes()); // Object terminator
    }

    /// Extracts all `0x0803` (`ChainAngle`) properties in degrees across the document tree.
    pub fn chain_angles_degrees(&self) -> Vec<f64> {
        let mut out = Vec::new();
        for obj in &self.root_objects {
            Self::collect_chain_angles(obj, &mut out);
        }
        out
    }

    fn collect_chain_angles(obj: &CdxObject, out: &mut Vec<f64>) {
        for ch in &obj.children {
            match ch {
                CdxChild::Property(p) => {
                    if let Some(deg) = p.as_chain_angle_degrees() {
                        out.push(deg);
                    }
                }
                CdxChild::Object(sub) => Self::collect_chain_angles(sub, out),
            }
        }
    }

    /// Extracts all chemical `0x8003` (`kCDXObj_Fragment`) objects into `(fragment_id, ChemMolecule)`.
    pub fn extract_molecules(&self) -> Vec<(u32, ChemMolecule)> {
        let mut out = Vec::new();
        for obj in &self.root_objects {
            Self::collect_cdx_fragments(obj, &mut out);
        }
        out
    }

    fn collect_cdx_fragments(obj: &CdxObject, out: &mut Vec<(u32, ChemMolecule)>) {
        if obj.tag == CDX_OBJ_FRAGMENT {
            let mut mol = ChemMolecule::new();
            let mut id_to_idx: HashMap<u32, usize> = HashMap::new();

            for ch in &obj.children {
                if let CdxChild::Object(node) = ch {
                    if node.tag == CDX_OBJ_NODE {
                        let mut node_type: i16 = 1; // Unspecified/Element
                        let mut z: i32 = 6; // Default Carbon
                        let mut has_element = false;
                        let mut charge = 0i32;
                        let mut isotope = 0i32;
                        let mut radical = 0i32;
                        let mut x = 0.0f64;
                        let mut y = 0.0f64;

                        for nch in &node.children {
                            if let CdxChild::Property(p) = nch {
                                match p.tag {
                                    CDX_PROP_NODE_TYPE => {
                                        node_type = p.as_i16().unwrap_or(1);
                                    }
                                    CDX_PROP_NODE_ELEMENT => {
                                        z = p.as_i16().unwrap_or(6) as i32;
                                        has_element = true;
                                    }
                                    CDX_PROP_ATOM_CHARGE => {
                                        charge = p.as_i8().unwrap_or(0) as i32;
                                    }
                                    CDX_PROP_ATOM_ISOTOPE => {
                                        isotope = p.as_i16().unwrap_or(0) as i32;
                                    }
                                    CDX_PROP_ATOM_RADICAL => {
                                        radical = p.as_i8().unwrap_or(0) as i32;
                                    }
                                    CDX_PROP_2D_POSITION => {
                                        if let Some((px, py)) = p.as_2d_position_xy_pt() {
                                            x = px;
                                            y = py;
                                        }
                                    }
                                    _ => {}
                                }
                            }
                        }

                        if !has_element && node_type != 1 {
                            z = 0;
                        }
                        let mut atom = ChemAtom::new(node.id, z);
                        atom.charge = charge;
                        atom.isotope = isotope;
                        atom.radical = radical;
                        atom.is_query = matches!(node_type, 2 | 3 | 6 | 7 | 8 | 9 | 11);
                        atom.x = x;
                        atom.y = y;
                        let idx = mol.add_atom(atom);
                        id_to_idx.insert(node.id, idx);
                    }
                }
            }

            for ch in &obj.children {
                if let CdxChild::Object(bond) = ch {
                    if bond.tag == CDX_OBJ_BOND {
                        let mut b_id: Option<u32> = None;
                        let mut e_id: Option<u32> = None;
                        let mut order_bits: u16 = 0x0001;
                        let mut display_code: u16 = 0;

                        for bch in &bond.children {
                            if let CdxChild::Property(p) = bch {
                                match p.tag {
                                    CDX_PROP_BOND_BEGIN => b_id = p.as_u32(),
                                    CDX_PROP_BOND_END => e_id = p.as_u32(),
                                    CDX_PROP_BOND_ORDER => {
                                        order_bits = p.as_u16().unwrap_or(0x0001);
                                    }
                                    CDX_PROP_BOND_DISPLAY => {
                                        display_code = p.as_u16().unwrap_or(0);
                                    }
                                    _ => {}
                                }
                            }
                        }

                        if let (Some(bi_id), Some(ei_id)) = (b_id, e_id) {
                            if let (Some(&bi), Some(&ei)) =
                                (id_to_idx.get(&bi_id), id_to_idx.get(&ei_id))
                            {
                                let order = match order_bits {
                                    0x0001 => 1,
                                    0x0002 => 2,
                                    0x0004 => 3,
                                    0x0008 => 4,
                                    0x0010 => 5, // Dative
                                    0x0040 => 7, // Hydrogen
                                    0x0080 => 6, // Delocalized / 1.5
                                    _ => 1,
                                };
                                let display = match display_code {
                                    0 => BondDisplay::Solid,
                                    1 => BondDisplay::Dash,
                                    2 => BondDisplay::Hash,
                                    3 => BondDisplay::WedgedHashBegin,
                                    4 => BondDisplay::WedgedHashEnd,
                                    5 => BondDisplay::Bold,
                                    6 => BondDisplay::WedgeBegin,
                                    7 => BondDisplay::WedgeEnd,
                                    8 => BondDisplay::Wavy,
                                    9 => BondDisplay::HollowWedgeBegin,
                                    10 => BondDisplay::HollowWedgeEnd,
                                    12 => BondDisplay::Dative,
                                    _ => BondDisplay::Solid,
                                };
                                mol.add_bond(bond.id, bi, ei, order, display);
                            }
                        }
                    }
                }
            }

            out.push((obj.id, mol));
        }

        for ch in &obj.children {
            if let CdxChild::Object(sub) = ch {
                Self::collect_cdx_fragments(sub, out);
            }
        }
    }

    /// Converts this `CdxDocument` into a `CdxmlDocument` DOM (pages, fragments, nodes, bonds, texts, settings).
    pub fn to_cdxml(&self) -> CdxmlDocument {
        let mut elements = Vec::new();
        let root_cdx = self.root_objects.first();
        let mut doc_attrs = Vec::new();
        if let Some(doc_obj) = root_cdx {
            for ch in &doc_obj.children {
                if let CdxChild::Property(p) = ch {
                    match p.tag {
                        CDX_PROP_CHAIN_ANGLE => {
                            if let Some(a) = p.as_chain_angle_degrees() {
                                doc_attrs.push(("ChainAngle".to_string(), format!("{}", a)));
                            }
                        }
                        CDX_PROP_BOND_SPACING => {
                            if let Some(sp) = p.as_u16() {
                                doc_attrs.push(("BondSpacing".to_string(), format!("{:.1}", (sp as f64) / 10.0)));
                            }
                        }
                        CDX_PROP_BOND_LENGTH => {
                            if let Some(l) = p.as_pt_1_65536() {
                                doc_attrs.push(("BondLength".to_string(), format!("{:.2}", l)));
                            }
                        }
                        CDX_PROP_BOLD_WIDTH => {
                            if let Some(w) = p.as_pt_1_65536() {
                                doc_attrs.push(("BoldWidth".to_string(), format!("{:.2}", w)));
                            }
                        }
                        CDX_PROP_LINE_WIDTH => {
                            if let Some(w) = p.as_pt_1_65536() {
                                doc_attrs.push(("LineWidth".to_string(), format!("{:.2}", w)));
                            }
                        }
                        CDX_PROP_MARGIN_WIDTH => {
                            if let Some(w) = p.as_pt_1_65536() {
                                doc_attrs.push(("MarginWidth".to_string(), format!("{:.2}", w)));
                            }
                        }
                        CDX_PROP_HASH_SPACING => {
                            if let Some(w) = p.as_pt_1_65536() {
                                doc_attrs.push(("HashSpacing".to_string(), format!("{:.2}", w)));
                            }
                        }
                        CDX_PROP_LABEL_STYLE => {
                            if p.payload.len() >= 6 {
                                let sz = u16::from_le_bytes([p.payload[4], p.payload[5]]) as f64 / 20.0;
                                if sz > 1.0 {
                                    doc_attrs.push(("LabelSize".to_string(), format!("{:.1}", sz)));
                                }
                            }
                        }
                        CDX_PROP_CAPTION_STYLE => {
                            if p.payload.len() >= 6 {
                                let sz = u16::from_le_bytes([p.payload[4], p.payload[5]]) as f64 / 20.0;
                                if sz > 1.0 {
                                    doc_attrs.push(("CaptionSize".to_string(), format!("{:.1}", sz)));
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
        elements.push(CdxmlElement {
            occurrence_id: 0,
            parent: None,
            tag: "CDXML".to_string(),
            attributes: doc_attrs,
            children: Vec::new(),
            self_closing: false,
        });

        if let Some(doc_obj) = root_cdx {
            for ch in &doc_obj.children {
                if let CdxChild::Object(sub) = ch {
                    Self::emit_cdxml_object(sub, 0, &mut elements);
                }
            }
        }

        let mut object_ids: HashMap<String, Vec<OccurrenceId>> = HashMap::new();
        for el in &elements {
            if let Some(id_val) = el.attr("id") {
                object_ids
                    .entry(id_val.to_string())
                    .or_default()
                    .push(el.occurrence_id);
            }
        }

        CdxmlDocument {
            prolog: "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n".to_string(),
            elements,
            root: 0,
            font_ids: HashMap::new(),
            object_ids,
        }
    }

    fn emit_cdxml_object(
        obj: &CdxObject,
        parent: OccurrenceId,
        elements: &mut Vec<CdxmlElement>,
    ) {
        let tag = match obj.tag {
            CDX_OBJ_PAGE => "page",
            CDX_OBJ_GROUP => "group",
            CDX_OBJ_FRAGMENT => "fragment",
            CDX_OBJ_NODE => "n",
            CDX_OBJ_BOND => "b",
            CDX_OBJ_TEXT => "t",
            CDX_OBJ_GRAPHIC => "graphic",
            CDX_OBJ_BRACKETED_GROUP => "bracketedgroup",
            _ => return,
        };

        let mut attrs = vec![("id".to_string(), obj.id.to_string())];
        let mut text_body: Option<String> = None;

        for ch in &obj.children {
            if let CdxChild::Property(p) = ch {
                match p.tag {
                    CDX_PROP_2D_POSITION => {
                        if let Some((x, y)) = p.as_2d_position_xy_pt() {
                            attrs.push(("p".to_string(), format!("{:.2} {:.2}", x, y)));
                        }
                    }
                    CDX_PROP_NODE_ELEMENT => {
                        if let Some(z) = p.as_i16() {
                            attrs.push(("Element".to_string(), z.to_string()));
                        }
                    }
                    CDX_PROP_ATOM_CHARGE => {
                        if let Some(c) = p.as_i8() {
                            attrs.push(("Charge".to_string(), c.to_string()));
                        }
                    }
                    CDX_PROP_ATOM_ISOTOPE => {
                        if let Some(iso) = p.as_i16() {
                            attrs.push(("Isotope".to_string(), iso.to_string()));
                        }
                    }
                    CDX_PROP_BOND_BEGIN => {
                        if let Some(b) = p.as_u32() {
                            attrs.push(("B".to_string(), b.to_string()));
                        }
                    }
                    CDX_PROP_BOND_END => {
                        if let Some(e) = p.as_u32() {
                            attrs.push(("E".to_string(), e.to_string()));
                        }
                    }
                    CDX_PROP_BOND_ORDER => {
                        if let Some(ord) = p.as_u16() {
                            let s = match ord {
                                0x0002 => "2",
                                0x0004 => "3",
                                0x0008 => "4",
                                0x0010 => "dative",
                                0x0040 => "hydrogen",
                                0x0080 => "1.5",
                                _ => "1",
                            };
                            attrs.push(("Order".to_string(), s.to_string()));
                        }
                    }
                    CDX_PROP_BOND_DISPLAY => {
                        if let Some(d) = p.as_u16() {
                            let s = match d {
                                1 => "Dash",
                                2 => "Hash",
                                3 => "WedgedHashBegin",
                                4 => "WedgedHashEnd",
                                5 => "Bold",
                                6 => "WedgeBegin",
                                7 => "WedgeEnd",
                                8 => "Wavy",
                                9 => "HollowWedgeBegin",
                                10 => "HollowWedgeEnd",
                                12 => "Dative",
                                _ => "Solid",
                            };
                            attrs.push(("Display".to_string(), s.to_string()));
                        }
                    }
                    CDX_PROP_TEXT => {
                        text_body = p.as_cdx_text();
                    }
                    _ => {}
                }
            }
        }

        let occ_id = elements.len();
        elements.push(CdxmlElement {
            occurrence_id: occ_id,
            parent: Some(parent),
            tag: tag.to_string(),
            attributes: attrs,
            children: Vec::new(),
            self_closing: false,
        });
        elements[parent]
            .children
            .push(CdxmlNodeChild::Element(occ_id));

        if let Some(txt) = text_body {
            let s_occ = elements.len();
            elements.push(CdxmlElement {
                occurrence_id: s_occ,
                parent: Some(occ_id),
                tag: "s".to_string(),
                attributes: Vec::new(),
                children: vec![CdxmlNodeChild::Text(txt)],
                self_closing: false,
            });
            elements[occ_id]
                .children
                .push(CdxmlNodeChild::Element(s_occ));
        }

        for ch in &obj.children {
            if let CdxChild::Object(sub) = ch {
                Self::emit_cdxml_object(sub, occ_id, elements);
            }
        }
        if elements[occ_id].children.is_empty() {
            elements[occ_id].self_closing = true;
        }
    }
}

// ============================================================================
// 2. CDXML Lossless DOM, Semantic Reference Resolver, & Roundtrip Serializer
// ============================================================================

/// Unique internal occurrence handle (`0..num_elements`), distinct from lexical `id` attribute.
pub type OccurrenceId = usize;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdxmlNodeChild {
    Element(OccurrenceId),
    Text(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdxmlElement {
    pub occurrence_id: OccurrenceId,
    pub parent: Option<OccurrenceId>,
    pub tag: String,
    /// Ordered `(name, value)` attributes preserving source order and unknown attributes.
    pub attributes: Vec<(String, String)>,
    pub children: Vec<CdxmlNodeChild>,
    pub self_closing: bool,
}

impl CdxmlElement {
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceResolution {
    Resolved(OccurrenceId),
    Ambiguous(Vec<OccurrenceId>),
    Dangling,
}

#[derive(Debug, Clone)]
pub struct CdxmlReferenceReport {
    pub duplicate_object_id_groups: BTreeMap<String, Vec<OccurrenceId>>,
    pub font_object_id_collisions: Vec<String>,
    pub bond_endpoints_total: usize,
    pub bond_endpoints_typed_local_resolved: usize,
    pub bond_endpoints_typed_local_ambiguous: usize,
    pub bond_endpoints_typed_local_dangling: usize,
    pub non_bond_references_resolved: usize,
    pub non_bond_references_unresolved: usize,
}

#[derive(Debug, Clone)]
pub struct CdxmlDocument {
    pub prolog: String,
    pub elements: Vec<CdxmlElement>,
    pub root: OccurrenceId,
    /// Separate namespace for `<font id="..."/>` in `<fonttable>` (74 corpus files collide with object IDs).
    pub font_ids: HashMap<String, Vec<OccurrenceId>>,
    /// Multi-valued object ID lookup (8 corpus files contain 34 duplicate non-font ID groups).
    pub object_ids: HashMap<String, Vec<OccurrenceId>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdxmlParserBudget {
    pub max_bytes: usize,
    pub max_elements: usize,
    pub max_depth: usize,
}

impl Default for CdxmlParserBudget {
    fn default() -> Self {
        Self {
            max_bytes: 64 * 1024 * 1024,
            max_elements: 2_000_000,
            max_depth: 256,
        }
    }
}

impl CdxmlDocument {
    pub fn parse(xml: &str) -> Result<Self, String> {
        Self::parse_with_budget(xml, CdxmlParserBudget::default())
    }

    pub fn parse_with_budget(xml: &str, budget: CdxmlParserBudget) -> Result<Self, String> {
        if xml.len() > budget.max_bytes {
            return Err("CDXML input exceeds max_bytes budget".to_string());
        }
        // Security: reject external entity declarations (<!ENTITY) while preserving <!DOCTYPE CDXML ...>
        if xml.contains("<!ENTITY") {
            return Err("External entity declarations (<!ENTITY) are prohibited".to_string());
        }
        let bytes = xml.as_bytes();
        let mut pos = 0usize;
        let len = bytes.len();

        // Capture prolog (<?xml ...?>, <!DOCTYPE ...>, whitespace/comments before root <CDXML)
        while pos < len {
            if bytes[pos] == b'<' {
                if pos + 2 <= len && (bytes[pos + 1] == b'?' || bytes[pos + 1] == b'!') {
                    if bytes[pos..].starts_with(b"<!--") {
                        if let Some(end) = xml[pos + 4..].find("-->") {
                            pos += 4 + end + 3;
                            continue;
                        } else {
                            return Err("Unterminated XML comment".to_string());
                        }
                    }
                    if let Some(end) = xml[pos..].find('>') {
                        pos += end + 1;
                        continue;
                    } else {
                        return Err("Unterminated XML prolog declaration".to_string());
                    }
                } else {
                    break;
                }
            } else {
                pos += 1;
            }
        }

        let prolog = xml[..pos].to_string();
        let mut elements: Vec<CdxmlElement> = Vec::new();
        let mut stack: Vec<OccurrenceId> = Vec::new();
        let mut root_id: Option<OccurrenceId> = None;

        while pos < len {
            if bytes[pos] == b'<' {
                if pos + 4 <= len && &bytes[pos..pos + 4] == b"<!--" {
                    if let Some(end) = xml[pos + 4..].find("-->") {
                        pos += 4 + end + 3;
                        continue;
                    } else {
                        return Err("Unterminated XML comment".to_string());
                    }
                }
                if pos + 2 <= len && bytes[pos + 1] == b'/' {
                    let close_end = xml[pos + 2..]
                        .find('>')
                        .ok_or_else(|| "Unterminated closing tag".to_string())?;
                    let close_tag = xml[pos + 2..pos + 2 + close_end].trim();
                    pos += 2 + close_end + 1;
                    let open_id = stack
                        .pop()
                        .ok_or_else(|| format!("Unexpected closing tag </{}>", close_tag))?;
                    if elements[open_id].tag != close_tag {
                        return Err(format!(
                            "Mismatched closing tag: expected </{}>, got </{}>",
                            elements[open_id].tag, close_tag
                        ));
                    }
                    if stack.is_empty() {
                        break;
                    }
                    continue;
                }

                // Opening or self-closing tag
                pos += 1; // skip '<'
                let tag_start = pos;
                while pos < len
                    && !bytes[pos].is_ascii_whitespace()
                    && bytes[pos] != b'/'
                    && bytes[pos] != b'>'
                {
                    pos += 1;
                }
                let tag_name = xml[tag_start..pos].to_string();
                if tag_name.is_empty() {
                    return Err(format!("Empty tag name at byte {}", tag_start));
                }

                let mut attrs = Vec::new();
                let mut self_closing = false;
                loop {
                    while pos < len && bytes[pos].is_ascii_whitespace() {
                        pos += 1;
                    }
                    if pos >= len {
                        return Err(format!("Unterminated tag <{}", tag_name));
                    }
                    if bytes[pos] == b'>' {
                        pos += 1;
                        break;
                    }
                    if pos + 2 <= len && &bytes[pos..pos + 2] == b"/>" {
                        self_closing = true;
                        pos += 2;
                        break;
                    }
                    let attr_start = pos;
                    while pos < len
                        && bytes[pos] != b'='
                        && !bytes[pos].is_ascii_whitespace()
                        && bytes[pos] != b'>'
                    {
                        pos += 1;
                    }
                    let attr_name = xml[attr_start..pos].to_string();
                    while pos < len && bytes[pos].is_ascii_whitespace() {
                        pos += 1;
                    }
                    if pos >= len || bytes[pos] != b'=' {
                        return Err(format!("Expected '=' after attribute {} in <{}>", attr_name, tag_name));
                    }
                    pos += 1;
                    while pos < len && bytes[pos].is_ascii_whitespace() {
                        pos += 1;
                    }
                    if pos >= len || (bytes[pos] != b'"' && bytes[pos] != b'\'') {
                        return Err(format!("Expected quoted value for attribute {}", attr_name));
                    }
                    let quote = bytes[pos];
                    pos += 1;
                    let val_start = pos;
                    while pos < len && bytes[pos] != quote {
                        pos += 1;
                    }
                    let raw_val = &xml[val_start..pos];
                    pos += 1;
                    attrs.push((attr_name, unescape_xml(raw_val)));
                }

                if elements.len() >= budget.max_elements {
                    return Err("CDXML element count exceeds max_elements budget".to_string());
                }
                if stack.len() >= budget.max_depth {
                    return Err("CDXML nesting exceeds max_depth budget".to_string());
                }
                let occ_id = elements.len();
                let parent = stack.last().copied();
                elements.push(CdxmlElement {
                    occurrence_id: occ_id,
                    parent,
                    tag: tag_name,
                    attributes: attrs,
                    children: Vec::new(),
                    self_closing,
                });
                if let Some(p) = parent {
                    elements[p].children.push(CdxmlNodeChild::Element(occ_id));
                } else {
                    root_id = Some(occ_id);
                }
                if !self_closing {
                    stack.push(occ_id);
                } else if stack.is_empty() {
                    break;
                }
            } else {
                let text_start = pos;
                while pos < len && bytes[pos] != b'<' {
                    pos += 1;
                }
                if let Some(&curr) = stack.last() {
                    let raw_txt = &xml[text_start..pos];
                    // Retain text inside `<s ...>` or `<t ...>` or non-whitespace text nodes
                    if elements[curr].tag == "s" || !raw_txt.trim().is_empty() {
                        elements[curr]
                            .children
                            .push(CdxmlNodeChild::Text(unescape_xml(raw_txt)));
                    }
                }
            }
        }

        let root = root_id.ok_or_else(|| "No root element found in CDXML".to_string())?;
        let mut font_ids: HashMap<String, Vec<OccurrenceId>> = HashMap::new();
        let mut object_ids: HashMap<String, Vec<OccurrenceId>> = HashMap::new();
        for el in &elements {
            if let Some(id_val) = el.attr("id") {
                if el.tag == "font" {
                    font_ids.entry(id_val.to_string()).or_default().push(el.occurrence_id);
                } else {
                    object_ids
                        .entry(id_val.to_string())
                        .or_default()
                        .push(el.occurrence_id);
                }
            }
        }

        Ok(Self {
            prolog,
            elements,
            root,
            font_ids,
            object_ids,
        })
    }

    /// Finds the nearest ancestor of `occ` whose tag equals `ancestor_tag`.
    pub fn ancestor(&self, mut occ: OccurrenceId, ancestor_tag: &str) -> Option<OccurrenceId> {
        while let Some(p) = self.elements[occ].parent {
            if self.elements[p].tag == ancestor_tag {
                return Some(p);
            }
            occ = p;
        }
        None
    }

    /// Resolves a bond endpoint (`B` or `E`) on bond `bond_occ` by restricting
    /// candidates to `<n>` elements within the same nearest `<fragment>` ancestor.
    /// Verified to resolve all 27,234 bond endpoints across all 655 shipped CDXML files uniquely.
    pub fn resolve_bond_endpoint(&self, bond_occ: OccurrenceId, target_id: &str) -> ReferenceResolution {
        let Some(candidates) = self.object_ids.get(target_id) else {
            return ReferenceResolution::Dangling;
        };
        let bond_frag = self.ancestor(bond_occ, "fragment");
        let typed_local: Vec<OccurrenceId> = candidates
            .iter()
            .copied()
            .filter(|&c| self.elements[c].tag == "n" && self.ancestor(c, "fragment") == bond_frag)
            .collect();
        match typed_local.len() {
            1 => ReferenceResolution::Resolved(typed_local[0]),
            0 => ReferenceResolution::Dangling,
            _ => ReferenceResolution::Ambiguous(typed_local),
        }
    }

    /// Runs the full reference & collision audit for parsed CDXML object IDs and cross-references.
    pub fn analyze_references(&self) -> CdxmlReferenceReport {
        let mut duplicate_object_id_groups = BTreeMap::new();
        for (id, occs) in &self.object_ids {
            if occs.len() > 1 {
                duplicate_object_id_groups.insert(id.clone(), occs.clone());
            }
        }
        let mut font_object_id_collisions: Vec<String> = self
            .font_ids
            .keys()
            .filter(|k| self.object_ids.contains_key(*k))
            .cloned()
            .collect();
        font_object_id_collisions.sort();

        let mut bond_endpoints_total = 0;
        let mut bond_endpoints_typed_local_resolved = 0;
        let mut bond_endpoints_typed_local_ambiguous = 0;
        let mut bond_endpoints_typed_local_dangling = 0;
        let mut non_bond_references_resolved = 0;
        let mut non_bond_references_unresolved = 0;

        const STEP_ATTRS: &[&str] = &[
            "ReactionStepReactants",
            "ReactionStepProducts",
            "ReactionStepArrows",
            "ReactionStepObjectsAboveArrow",
            "ReactionStepObjectsBelowArrow",
        ];

        for el in &self.elements {
            match el.tag.as_str() {
                "b" => {
                    for attr in ["B", "E"] {
                        if let Some(val) = el.attr(attr) {
                            for token in val.split_whitespace() {
                                bond_endpoints_total += 1;
                                match self.resolve_bond_endpoint(el.occurrence_id, token) {
                                    ReferenceResolution::Resolved(_) => {
                                        bond_endpoints_typed_local_resolved += 1
                                    }
                                    ReferenceResolution::Ambiguous(_) => {
                                        bond_endpoints_typed_local_ambiguous += 1
                                    }
                                    ReferenceResolution::Dangling => {
                                        bond_endpoints_typed_local_dangling += 1
                                    }
                                }
                            }
                        }
                    }
                }
                "step" => {
                    for &attr in STEP_ATTRS {
                        if let Some(val) = el.attr(attr) {
                            for token in val.split_whitespace() {
                                if self.object_ids.get(token).map_or(0, |v| v.len()) == 1 {
                                    non_bond_references_resolved += 1;
                                } else {
                                    non_bond_references_unresolved += 1;
                                }
                            }
                        }
                    }
                }
                "bracketedgroup" => {
                    if let Some(val) = el.attr("BracketedObjectIDs") {
                        for token in val.split_whitespace() {
                            if self.object_ids.get(token).map_or(0, |v| v.len()) == 1 {
                                non_bond_references_resolved += 1;
                            } else {
                                non_bond_references_unresolved += 1;
                            }
                        }
                    }
                }
                "represent" => {
                    if let Some(val) = el.attr("object") {
                        for token in val.split_whitespace() {
                            if self.object_ids.get(token).map_or(0, |v| v.len()) == 1 {
                                non_bond_references_resolved += 1;
                            } else {
                                non_bond_references_unresolved += 1;
                            }
                        }
                    }
                }
                "graphic" => {
                    if let Some(val) = el.attr("SupersededBy") {
                        for token in val.split_whitespace() {
                            if self.object_ids.get(token).map_or(0, |v| v.len()) == 1 {
                                non_bond_references_resolved += 1;
                            } else {
                                non_bond_references_unresolved += 1;
                            }
                        }
                    }
                }
                _ => {}
            }
        }

        CdxmlReferenceReport {
            duplicate_object_id_groups,
            font_object_id_collisions,
            bond_endpoints_total,
            bond_endpoints_typed_local_resolved,
            bond_endpoints_typed_local_ambiguous,
            bond_endpoints_typed_local_dangling,
            non_bond_references_resolved,
            non_bond_references_unresolved,
        }
    }

    /// Extracts all chemical `<fragment>` objects into `ChemMolecule` instances for `chiralcraft-chem`.
    pub fn extract_molecules(&self) -> Vec<(OccurrenceId, ChemMolecule)> {
        let mut out = Vec::new();
        for el in &self.elements {
            if el.tag == "fragment" {
                let mut mol = ChemMolecule::new();
                let mut occ_to_atom_idx: HashMap<OccurrenceId, usize> = HashMap::new();

                for ch in &el.children {
                    if let CdxmlNodeChild::Element(child_occ) = *ch {
                        let node = &self.elements[child_occ];
                        if node.tag == "n" {
                            let id = node
                                .attr("id")
                                .and_then(|s| s.parse::<u32>().ok())
                                .unwrap_or(child_occ as u32);
                            let node_type = node.attr("NodeType");
                            let z = node
                                .attr("Element")
                                .and_then(|s| s.parse::<i32>().ok())
                                .unwrap_or_else(|| {
                                    if node_type.is_some() {
                                        0
                                    } else {
                                        6 // Default unlabeled node in CDXML is Carbon (Z=6)
                                    }
                                });
                            let mut atom = ChemAtom::new(id, z);
                            if matches!(
                                node_type,
                                Some("GenericNickname")
                                    | Some("ElementList")
                                    | Some("AnonymousAlternativeGroup")
                                    | Some("MultiAttachment")
                                    | Some("VariableAttachment")
                                    | Some("ExternalConnectionPoint")
                            ) {
                                atom.is_query = true;
                            } else if node_type == Some("Unspecified") {
                                atom.is_unspecified_label = true;
                            }
                            atom.visible = node.attr("Visible") != Some("no");
                            atom.ignore_warnings = node.attr("IgnoreWarnings") == Some("yes");
                            if let Some(chg) = node.attr("Charge").and_then(|s| s.parse::<i32>().ok()) {
                                atom.charge = chg;
                            }
                            if let Some(iso) = node.attr("Isotope").and_then(|s| s.parse::<i32>().ok()) {
                                atom.isotope = iso;
                            }
                            if let Some(rad_str) = node.attr("Radical") {
                                atom.radical = match rad_str {
                                    "Doublet" | "1" => 1,
                                    "Singlet" | "2" => 2,
                                    "Triplet" | "3" => 3,
                                    _ => rad_str.parse::<i32>().unwrap_or(0),
                                };
                            }
                            if let Some(val_ovr) = node.attr("Valence").and_then(|s| s.parse::<i32>().ok()) {
                                atom.valence_override = Some(val_ovr);
                            }
                            if let Some(nh) = node.attr("NumHydrogens").and_then(|s| s.parse::<i32>().ok()) {
                                atom.explicit_num_hydrogens = Some(nh);
                            }
                            if let Some(p_str) = node.attr("p") {
                                let mut parts = p_str.split_whitespace();
                                if let (Some(xs), Some(ys)) = (parts.next(), parts.next()) {
                                    if let (Ok(x), Ok(y)) = (xs.parse::<f64>(), ys.parse::<f64>()) {
                                        atom.x = x;
                                        atom.y = y;
                                    }
                                }
                            }
                            if let Some(xyz_str) = node.attr("xyz") {
                                let mut parts = xyz_str.split_whitespace();
                                if let (Some(xs), Some(ys), Some(zs)) =
                                    (parts.next(), parts.next(), parts.next())
                                {
                                    if let (Ok(x), Ok(y), Ok(zc)) =
                                        (xs.parse::<f64>(), ys.parse::<f64>(), zs.parse::<f64>())
                                    {
                                        atom.x = x;
                                        atom.y = y;
                                        atom.z_coord = zc;
                                        mol.has_3d_coords = true;
                                    }
                                }
                            }
                            if let Some(clbl) = node
                                .attr("CustomLabel")
                                .or_else(|| node.attr("GenericNickname"))
                            {
                                if !clbl.is_empty() {
                                    atom.label = Some(clbl.to_string());
                                }
                            }
                            let idx = mol.add_atom(atom);
                            occ_to_atom_idx.insert(child_occ, idx);
                        }
                    }
                }

                for ch in &el.children {
                    if let CdxmlNodeChild::Element(child_occ) = *ch {
                        let bond = &self.elements[child_occ];
                        if bond.tag == "b" {
                            let id = bond
                                .attr("id")
                                .and_then(|s| s.parse::<u32>().ok())
                                .unwrap_or(child_occ as u32);
                            let b_ref = bond.attr("B").unwrap_or("");
                            let e_ref = bond.attr("E").unwrap_or("");
                            if let (
                                ReferenceResolution::Resolved(b_occ),
                                ReferenceResolution::Resolved(e_occ),
                            ) = (
                                self.resolve_bond_endpoint(child_occ, b_ref),
                                self.resolve_bond_endpoint(child_occ, e_ref),
                            ) {
                                if let (Some(&bi), Some(&ei)) =
                                    (occ_to_atom_idx.get(&b_occ), occ_to_atom_idx.get(&e_occ))
                                {
                                    let order = match bond.attr("Order") {
                                        None | Some("1") => 1,
                                        Some("2") => 2,
                                        Some("3") => 3,
                                        Some("4") => 4,
                                        Some("1.5") | Some("2.5") | Some("aromatic") => 6,
                                        Some("dative") => 5,
                                        Some("hydrogen") => 7,
                                        _ => 1,
                                    };
                                    let display = match bond.attr("Display") {
                                        None | Some("Solid") => BondDisplay::Solid,
                                        Some("Dash") => BondDisplay::Dash,
                                        Some("Hash") => BondDisplay::Hash,
                                        Some("WedgedHashBegin") => BondDisplay::WedgedHashBegin,
                                        Some("WedgedHashEnd") => BondDisplay::WedgedHashEnd,
                                        Some("Bold") => BondDisplay::Bold,
                                        Some("WedgeBegin") => BondDisplay::WedgeBegin,
                                        Some("WedgeEnd") => BondDisplay::WedgeEnd,
                                        Some("Wavy") => BondDisplay::Wavy,
                                        Some("HollowWedgeBegin") => BondDisplay::HollowWedgeBegin,
                                        Some("HollowWedgeEnd") => BondDisplay::HollowWedgeEnd,
                                        Some("Dative") => BondDisplay::Dative,
                                        _ => BondDisplay::Solid,
                                    };
                                    let b_idx = mol.add_bond(id, bi, ei, order, display);
                                    mol.bonds[b_idx].visible = bond.attr("Visible") != Some("no");
                                    mol.bonds[b_idx].ignore_warnings =
                                        bond.attr("IgnoreWarnings") == Some("yes");
                                }
                            }
                        }
                    }
                }

                out.push((el.occurrence_id, mol));
            }
        }
        out
    }

    /// Serializes the `CdxmlDocument` back to XML string, preserving all attributes,
    /// absent vs explicit `Order`, duplicate IDs, and unknown elements.
    pub fn to_xml_string(&self) -> String {
        let mut out = String::with_capacity(4096);
        out.push_str(&self.prolog);
        self.write_element(self.root, &mut out);
        out
    }

    /// Converts this `CdxmlDocument` into a standard 22-byte-header binary `CdxDocument` (`VjCD0100`).
    pub fn to_cdx(&self) -> CdxDocument {
        let mut header_bytes = Vec::with_capacity(22);
        header_bytes.extend_from_slice(CDX_MAGIC_PREFIX);
        header_bytes.extend_from_slice(&[0u8; 10]);

        let mut total_objects = 0usize;
        let mut total_properties = 0usize;
        let mut total_object_terminators = 0usize;

        fn prop(tag: u16, payload: Vec<u8>, props: &mut usize) -> CdxChild {
            *props += 1;
            CdxChild::Property(CdxProperty {
                offset: 0,
                tag,
                extended_length: payload.len() >= 0xffff,
                payload,
            })
        }

        fn gather_text(doc: &CdxmlDocument, occ: OccurrenceId, out: &mut String) {
            if let Some(el) = doc.elements.get(occ) {
                for ch in &el.children {
                    match ch {
                        CdxmlNodeChild::Text(t) => out.push_str(t),
                        CdxmlNodeChild::Element(cid) => gather_text(doc, *cid, out),
                    }
                }
            }
        }

        fn build_obj(
            doc: &CdxmlDocument,
            occ: OccurrenceId,
            objs: &mut usize,
            props: &mut usize,
            terms: &mut usize,
        ) -> Option<CdxObject> {
            let el = doc.elements.get(occ)?;
            let tag = match el.tag.as_str() {
                "CDXML" => CDX_OBJ_DOCUMENT,
                "page" => CDX_OBJ_PAGE,
                "group" => CDX_OBJ_GROUP,
                "fragment" => CDX_OBJ_FRAGMENT,
                "n" => CDX_OBJ_NODE,
                "b" => CDX_OBJ_BOND,
                "t" => CDX_OBJ_TEXT,
                "graphic" => CDX_OBJ_GRAPHIC,
                "bracketedgroup" => CDX_OBJ_BRACKETED_GROUP,
                _ => return None,
            };
            let id = el
                .attr("id")
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            let mut children = Vec::new();

            if tag == CDX_OBJ_DOCUMENT {
                if let Some(bl) = el.attr("BondLength").and_then(|s| s.parse::<f64>().ok()) {
                    let v = (bl * 65536.0).round() as i32;
                    children.push(prop(CDX_PROP_BOND_LENGTH, v.to_le_bytes().to_vec(), props));
                }
                if let Some(lw) = el.attr("LineWidth").and_then(|s| s.parse::<f64>().ok()) {
                    let v = (lw * 65536.0).round() as i32;
                    children.push(prop(CDX_PROP_LINE_WIDTH, v.to_le_bytes().to_vec(), props));
                }
                if let Some(ca) = el.attr("ChainAngle").and_then(|s| s.parse::<f64>().ok()) {
                    let v = (ca * 65536.0).round() as i32;
                    children.push(prop(CDX_PROP_CHAIN_ANGLE, v.to_le_bytes().to_vec(), props));
                }
            }

            if let Some(p_str) = el.attr("p") {
                let nums: Vec<f64> = p_str
                    .split_whitespace()
                    .filter_map(|t| t.parse().ok())
                    .collect();
                if nums.len() >= 2 {
                    // CDX 2D Position stores Y before X in 1/65536 pt!
                    let y = (nums[1] * 65536.0).round() as i32;
                    let x = (nums[0] * 65536.0).round() as i32;
                    let mut pl = Vec::with_capacity(8);
                    pl.extend_from_slice(&y.to_le_bytes());
                    pl.extend_from_slice(&x.to_le_bytes());
                    children.push(prop(CDX_PROP_2D_POSITION, pl, props));
                }
            }
            if let Some(z) = el.attr("Element").and_then(|s| s.parse::<i16>().ok()) {
                children.push(prop(CDX_PROP_NODE_ELEMENT, z.to_le_bytes().to_vec(), props));
            }
            if let Some(c) = el.attr("Charge").and_then(|s| s.parse::<i8>().ok()) {
                children.push(prop(CDX_PROP_ATOM_CHARGE, vec![c as u8], props));
            }
            if let Some(iso) = el.attr("Isotope").and_then(|s| s.parse::<i16>().ok()) {
                children.push(prop(CDX_PROP_ATOM_ISOTOPE, iso.to_le_bytes().to_vec(), props));
            }
            if let Some(b) = el.attr("B").and_then(|s| s.parse::<u32>().ok()) {
                children.push(prop(CDX_PROP_BOND_BEGIN, b.to_le_bytes().to_vec(), props));
            }
            if let Some(e) = el.attr("E").and_then(|s| s.parse::<u32>().ok()) {
                children.push(prop(CDX_PROP_BOND_END, e.to_le_bytes().to_vec(), props));
            }
            if let Some(ord) = el.attr("Order") {
                let code: u16 = match ord {
                    "2" => 0x0002,
                    "3" => 0x0004,
                    "4" => 0x0008,
                    "dative" => 0x0010,
                    "hydrogen" => 0x0040,
                    "1.5" => 0x0080,
                    _ => 0x0001,
                };
                children.push(prop(CDX_PROP_BOND_ORDER, code.to_le_bytes().to_vec(), props));
            }
            if let Some(disp) = el.attr("Display") {
                let code: u16 = match disp {
                    "Dash" => 1,
                    "Hash" => 2,
                    "WedgedHashBegin" => 3,
                    "WedgedHashEnd" => 4,
                    "Bold" => 5,
                    "WedgeBegin" => 6,
                    "WedgeEnd" => 7,
                    "Wavy" => 8,
                    "HollowWedgeBegin" => 9,
                    "HollowWedgeEnd" => 10,
                    "Dative" => 12,
                    _ => 0,
                };
                children.push(prop(CDX_PROP_BOND_DISPLAY, code.to_le_bytes().to_vec(), props));
            }
            if tag == CDX_OBJ_TEXT {
                let mut txt = String::new();
                gather_text(doc, occ, &mut txt);
                if !txt.is_empty() {
                    let mut pl = vec![0u8, 0u8]; // 0 style runs
                    pl.extend_from_slice(txt.as_bytes());
                    children.push(prop(CDX_PROP_TEXT, pl, props));
                }
            }

            for ch in &el.children {
                if let CdxmlNodeChild::Element(cid) = ch {
                    if let Some(sub) = build_obj(doc, *cid, objs, props, terms) {
                        children.push(CdxChild::Object(sub));
                    }
                }
            }

            *objs += 1;
            *terms += 1;
            Some(CdxObject {
                offset: 0,
                tag,
                id,
                children,
            })
        }

        let mut root_objects = Vec::new();
        if let Some(root_obj) = build_obj(
            self,
            self.root,
            &mut total_objects,
            &mut total_properties,
            &mut total_object_terminators,
        ) {
            root_objects.push(root_obj);
        }

        CdxDocument {
            header_bytes,
            root_objects,
            total_objects,
            total_properties,
            total_object_terminators,
        }
    }

    fn write_element(&self, occ: OccurrenceId, out: &mut String) {
        let el = &self.elements[occ];
        out.push('<');
        out.push_str(&el.tag);
        for (k, v) in &el.attributes {
            out.push(' ');
            out.push_str(k);
            out.push_str("=\"");
            out.push_str(&escape_xml_attr(v));
            out.push('"');
        }
        if el.self_closing && el.children.is_empty() {
            out.push_str("/>");
            return;
        }
        out.push('>');
        for ch in &el.children {
            match ch {
                CdxmlNodeChild::Element(c_occ) => self.write_element(*c_occ, out),
                CdxmlNodeChild::Text(t) => out.push_str(&escape_xml_text(t)),
            }
        }
        out.push_str("</");
        out.push_str(&el.tag);
        out.push('>');
    }
}

fn unescape_xml(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        rest = &rest[amp + 1..];
        if let Some(semi) = rest.find(';') {
            let ent = &rest[..semi];
            rest = &rest[semi + 1..];
            match ent {
                "amp" => out.push('&'),
                "lt" => out.push('<'),
                "gt" => out.push('>'),
                "quot" => out.push('"'),
                "apos" => out.push('\''),
                _ if ent.starts_with("#x") || ent.starts_with("#X") => {
                    if let Ok(code) = u32::from_str_radix(&ent[2..], 16) {
                        if let Some(ch) = char::from_u32(code) {
                            out.push(ch);
                        }
                    }
                }
                _ if ent.starts_with('#') => {
                    if let Ok(code) = ent[1..].parse::<u32>() {
                        if let Some(ch) = char::from_u32(code) {
                            out.push(ch);
                        }
                    }
                }
                _ => {
                    out.push('&');
                    out.push_str(ent);
                    out.push(';');
                }
            }
        } else {
            out.push('&');
        }
    }
    out.push_str(rest);
    out
}


fn escape_xml_attr(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}


fn escape_xml_text(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            _ => out.push(ch),
        }
    }
    out
}


/// Extracted reaction scheme & steps from `<scheme>` / `<step>` elements in CDXML.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedReactionStep {
    pub step_id: Option<String>,
    pub reactant_ids: Vec<String>,
    pub product_ids: Vec<String>,
    pub arrow_ids: Vec<String>,
    pub above_arrow_ids: Vec<String>,
    pub below_arrow_ids: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedReactionScheme {
    pub scheme_id: Option<String>,
    pub steps: Vec<ExtractedReactionStep>,
}

impl CdxmlDocument {
    /// Extracts all `<scheme>` and `<step>` reaction structures from the CDXML document.
    pub fn extract_reactions(&self) -> Vec<ExtractedReactionScheme> {
        let mut schemes = Vec::new();
        for el in &self.elements {
            if el.tag == "scheme" {
                let mut steps = Vec::new();
                for ch in &el.children {
                    if let CdxmlNodeChild::Element(child_occ) = *ch {
                        let st = &self.elements[child_occ];
                        if st.tag == "step" {
                            let split_ids = |attr: &str| -> Vec<String> {
                                st.attr(attr)
                                    .map(|v| v.split_whitespace().map(|s| s.to_string()).collect())
                                    .unwrap_or_default()
                            };
                            steps.push(ExtractedReactionStep {
                                step_id: st.attr("id").map(|s| s.to_string()),
                                reactant_ids: split_ids("ReactionStepReactants"),
                                product_ids: split_ids("ReactionStepProducts"),
                                arrow_ids: split_ids("ReactionStepArrows"),
                                above_arrow_ids: split_ids("ReactionStepObjectsAboveArrow"),
                                below_arrow_ids: split_ids("ReactionStepObjectsBelowArrow"),
                            });
                        }
                    }
                }
                schemes.push(ExtractedReactionScheme {
                    scheme_id: el.attr("id").map(|s| s.to_string()),
                    steps,
                });
            }
        }
        schemes
    }

    /// Extracts all `<bracketedgroup>` elements, attached bracket graphics, and crossing bonds.
    pub fn extract_bracketed_groups(&self) -> Vec<ExtractedBracketedGroup> {
        let mut groups = Vec::new();
        for el in &self.elements {
            if el.tag == "bracketedgroup" {
                let mut attachments = Vec::new();
                for ch in &el.children {
                    if let CdxmlNodeChild::Element(att_occ) = *ch {
                        let att = &self.elements[att_occ];
                        if att.tag == "bracketattachment" {
                            let mut crossing_bonds = Vec::new();
                            for bch in &att.children {
                                if let CdxmlNodeChild::Element(cb_occ) = *bch {
                                    let cb = &self.elements[cb_occ];
                                    if cb.tag == "crossingbond" {
                                        crossing_bonds.push(ExtractedCrossingBond {
                                            id: cb.attr("id").map(|s| s.to_string()),
                                            bond_id: cb.attr("BondID").map(|s| s.to_string()),
                                            inner_atom_id: cb.attr("InnerAtomID").map(|s| s.to_string()),
                                        });
                                    }
                                }
                            }
                            attachments.push(ExtractedBracketAttachment {
                                id: att.attr("id").map(|s| s.to_string()),
                                graphic_id: att.attr("GraphicID").map(|s| s.to_string()),
                                crossing_bonds,
                            });
                        }
                    }
                }
                let bracketed_object_ids = el
                    .attr("BracketedObjectIDs")
                    .map(|v| v.split_whitespace().map(|s| s.to_string()).collect())
                    .unwrap_or_default();
                groups.push(ExtractedBracketedGroup {
                    id: el.attr("id").map(|s| s.to_string()),
                    bracket_usage: el.attr("BracketUsage").map(|s| s.to_string()),
                    repeat_count: el.attr("RepeatCount").map(|s| s.to_string()),
                    flip_pattern: el.attr("FlipPattern").map(|s| s.to_string()),
                    bracketed_object_ids,
                    attachments,
                });
            }
        }
        groups
    }

    fn collect_text_rec(&self, occ: OccurrenceId, out: &mut String) {
        for ch in &self.elements[occ].children {
            match ch {
                CdxmlNodeChild::Text(t) => out.push_str(t),
                CdxmlNodeChild::Element(c_occ) => self.collect_text_rec(*c_occ, out),
            }
        }
    }

    /// Extracts `<chemicalproperty ChemicalPropertyType="23">` (`Object > Attach Data...`),
    /// resolving linked `BasisObjects`, display caption text (`D = 15.9`), and `<objecttag Name="attachedData">`.
    pub fn extract_attached_data(&self) -> Vec<ExtractedAttachedData> {
        let mut caption_by_id: HashMap<&str, String> = HashMap::new();
        let mut tag_by_cp_id: HashMap<&str, (Option<String>, String)> = HashMap::new();

        for (occ, el) in self.elements.iter().enumerate() {
            if el.tag == "t" {
                if let Some(id) = el.attr("id") {
                    let mut txt = String::new();
                    self.collect_text_rec(occ, &mut txt);
                    caption_by_id.insert(id, txt);
                }
            } else if el.tag == "objecttag" && el.attr("Name") == Some("attachedData") {
                if let Some(cp_id) = el.attr("Value") {
                    let mut tag_txt = String::new();
                    self.collect_text_rec(occ, &mut tag_txt);
                    let offset = el.attr("PositioningOffset").map(|s| s.to_string());
                    tag_by_cp_id.insert(cp_id, (offset, tag_txt));
                }
            }
        }

        let mut out = Vec::new();
        for el in &self.elements {
            if el.tag == "chemicalproperty" && el.attr("ChemicalPropertyType") == Some("23") {
                let id = el.attr("id").map(|s| s.to_string());
                let field_name = el.attr("Name").unwrap_or("").to_string();
                let display_id = el.attr("ChemicalPropertyDisplayID");
                let display_text = display_id
                    .and_then(|did| caption_by_id.get(did))
                    .cloned()
                    .unwrap_or_default();
                let basis_objects = el
                    .attr("BasisObjects")
                    .map(|v| v.split_whitespace().map(|s| s.to_string()).collect())
                    .unwrap_or_default();
                let (positioning_offset, tag) = id
                    .as_deref()
                    .and_then(|cid| tag_by_cp_id.get(cid))
                    .map(|(off, tg)| (off.clone(), Some(tg.clone()).filter(|s| !s.is_empty())))
                    .unwrap_or((None, None));
                out.push(ExtractedAttachedData {
                    id,
                    field_name,
                    tag,
                    display_text,
                    positioning_offset,
                    basis_objects,
                });
            }
        }
        out
    }

    /// Extracts `<annotation Keyword="..." Content="..."/>` elements (`Object > Annotate...`).
    pub fn extract_annotations(&self) -> Vec<ExtractedAnnotation> {
        self.elements
            .iter()
            .filter(|el| el.tag == "annotation")
            .map(|el| ExtractedAnnotation {
                parent_occ: el.parent,
                keyword: el.attr("Keyword").unwrap_or("").to_string(),
                content: el.attr("Content").unwrap_or("").to_string(),
            })
            .collect()
    }

    /// Extracts `<objecttag Name="stereo">` indicator captions (`Object > Show Stereochemistry`).
    pub fn extract_stereo_indicators(&self) -> Vec<ExtractedStereoIndicator> {
        let mut out = Vec::new();
        for (occ, el) in self.elements.iter().enumerate() {
            if el.tag == "objecttag" && el.attr("Name") == Some("stereo") {
                let mut text = String::new();
                self.collect_text_rec(occ, &mut text);
                let parent_id = el
                    .parent
                    .and_then(|p| self.elements.get(p))
                    .and_then(|pe| pe.attr("id"))
                    .map(str::to_string);
                out.push(ExtractedStereoIndicator {
                    parent_occ: el.parent,
                    parent_id,
                    text,
                });
            }
        }
        out
    }

    /// Extracts `<graphic GraphicType="Oval" OvalType="Circle">` aromatic ring circles (`Toggle Aromatic Display`).
    pub fn extract_aromatic_circles(&self) -> Vec<ExtractedAromaticCircle> {
        let mut out = Vec::new();
        for (occ, el) in self.elements.iter().enumerate() {
            if el.tag == "graphic"
                && el.attr("GraphicType") == Some("Oval")
                && el.attr("OvalType") == Some("Circle")
            {
                let mut cx = 0.0;
                let mut cy = 0.0;
                let mut radius = 0.0;
                if let Some(c_str) = el.attr("Center3D") {
                    let mut it = c_str.split_whitespace();
                    cx = it.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
                    cy = it.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
                    if let Some(m_str) = el.attr("MajorAxisEnd3D") {
                        let mut mit = m_str.split_whitespace();
                        let mx = mit.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(cx);
                        let my = mit.next().and_then(|s| s.parse::<f64>().ok()).unwrap_or(cy);
                        radius = ((mx - cx).powi(2) + (my - cy).powi(2)).sqrt();
                    }
                } else if let Some(bb_str) = el.attr("BoundingBox") {
                    let vals: Vec<f64> = bb_str
                        .split_whitespace()
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect();
                    if vals.len() == 4 {
                        cx = 0.5 * (vals[0] + vals[2]);
                        cy = 0.5 * (vals[1] + vals[3]);
                        radius = 0.5 * (vals[2] - vals[0]).abs();
                    }
                }
                out.push(ExtractedAromaticCircle {
                    fragment_occ: self.ancestor(occ, "fragment"),
                    id: el.attr("id").map(str::to_string),
                    center_x: cx,
                    center_y: cy,
                    radius,
                });
            }
        }
        out
    }

    /// Executes `Structure > Clear Reaction Map` in CDXML: removes `ReactionStepAtomMap` from
    /// `<step>` elements and `RxnParticipation` from `<b>` elements.
    pub fn clear_reaction_maps(&mut self) {
        for el in &mut self.elements {
            if el.tag == "step" {
                el.attributes.retain(|(k, _)| k != "ReactionStepAtomMap");
            } else if el.tag == "b" {
                el.attributes.retain(|(k, _)| k != "RxnParticipation");
            }
        }
    }

    /// Executes `Object > Flatten` in CDXML: removes `xyz` 3D coordinates from `<n>` nodes
    /// while preserving/ensuring 2D `p` coordinates.
    pub fn flatten_3d(&mut self) {
        for el in &mut self.elements {
            if el.tag == "n" {
                if let Some(xyz) = el.attr("xyz").map(str::to_string) {
                    if el.attr("p").is_none() {
                        let mut it = xyz.split_whitespace();
                        if let (Some(x), Some(y)) = (
                            it.next().and_then(|s| s.parse::<f64>().ok()),
                            it.next().and_then(|s| s.parse::<f64>().ok()),
                        ) {
                            el.attributes
                                .push(("p".to_string(), format!("{x:.2} {y:.2}")));
                        }
                    }
                    el.attributes.retain(|(k, _)| k != "xyz");
                }
            }
        }
    }
}

/// Extracted `<annotation Keyword="..." Content="..."/>` item from CDXML (`Object > Annotate...`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedAnnotation {
    pub parent_occ: Option<OccurrenceId>,
    pub keyword: String,
    pub content: String,
}

/// Extracted `<objecttag Name="stereo">` label from CDXML (`Object > Show Stereochemistry`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedStereoIndicator {
    pub parent_occ: Option<OccurrenceId>,
    pub parent_id: Option<String>,
    pub text: String,
}

/// Extracted `<graphic GraphicType="Oval" OvalType="Circle">` ring circle from CDXML (`Toggle Aromatic Display`).
#[derive(Debug, Clone, PartialEq)]
pub struct ExtractedAromaticCircle {
    pub fragment_occ: Option<OccurrenceId>,
    pub id: Option<String>,
    pub center_x: f64,
    pub center_y: f64,
    pub radius: f64,
}

/// Extracted `Object > Attach Data...` item (`ChemicalPropertyType="23"`) from CDXML (`C04`/`C12`/`C21`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedAttachedData {
    pub id: Option<String>,
    pub field_name: String,
    pub tag: Option<String>,
    pub display_text: String,
    pub positioning_offset: Option<String>,
    pub basis_objects: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedCrossingBond {
    pub id: Option<String>,
    pub bond_id: Option<String>,
    pub inner_atom_id: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedBracketAttachment {
    pub id: Option<String>,
    pub graphic_id: Option<String>,
    pub crossing_bonds: Vec<ExtractedCrossingBond>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtractedBracketedGroup {
    pub id: Option<String>,
    pub bracket_usage: Option<String>,
    pub repeat_count: Option<String>,
    pub flip_pattern: Option<String>,
    pub bracketed_object_ids: Vec<String>,
    pub attachments: Vec<ExtractedBracketAttachment>,
}

// ============================================================================
// 3. MDL Molfile V2000 Interchange (Import & Export)
// ============================================================================

#[inline]
fn mdl_rad_to_chem(mdl: i32) -> i32 {
    match mdl {
        2 => 1, // MDL Doublet (2) -> ChiralCraft Doublet (1)
        1 => 2, // MDL Singlet (1) -> ChiralCraft Singlet (2)
        3 => 3, // Triplet (3)
        _ => 0,
    }
}

#[inline]
fn chem_rad_to_mdl(chem: i32) -> i32 {
    match chem {
        1 => 2, // ChiralCraft Doublet (1) -> MDL Doublet (2)
        2 => 1, // ChiralCraft Singlet (2) -> MDL Singlet (1)
        3 => 3, // Triplet (3)
        _ => 0,
    }
}

/// Parses an MDL Molfile V2000 string into a `ChemMolecule`.
pub fn molfile_v2000_to_molecule(molfile: &str) -> Result<ChemMolecule, String> {
    let lines: Vec<&str> = molfile.lines().collect();
    if lines.len() < 4 {
        return Err("Molfile V2000 truncated header".to_string());
    }
    let counts = lines[3];
    if counts.len() < 6 {
        return Err("Molfile V2000 invalid counts line".to_string());
    }
    let n_atoms: usize = counts[0..3]
        .trim()
        .parse()
        .map_err(|_| "Invalid atom count in V2000 header".to_string())?;
    let n_bonds: usize = counts[3..6]
        .trim()
        .parse()
        .map_err(|_| "Invalid bond count in V2000 header".to_string())?;
    if lines.len() < 4 + n_atoms + n_bonds {
        return Err("Molfile V2000 truncated atom/bond block".to_string());
    }

    let mut mol = ChemMolecule::new();
    for i in 0..n_atoms {
        let line = lines[4 + i];
        let x: f64 = line.get(0..10).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
        let y: f64 = line.get(10..20).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
        let zc: f64 = line.get(20..30).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
        let sym = line.get(31..34).unwrap_or("C").trim();
        let chg_code: i32 = line.get(36..39).and_then(|s| s.trim().parse().ok()).unwrap_or(0);
        let z = if let Some(el) = get_element_by_symbol(sym) {
            el.z
        } else if sym == "D" {
            225
        } else if sym == "T" {
            226
        } else {
            0
        };
        let mut atom = ChemAtom::new(i as u32 + 1, z);
        atom.x = x;
        atom.y = -y; // Molfile Y-up -> ChiralCraft Y-down
        atom.z_coord = zc;
        if zc.abs() > 1e-6 {
            mol.has_3d_coords = true;
        }
        atom.charge = match chg_code {
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            _ => 0,
        };
        if chg_code == 4 {
            atom.radical = 1; // Doublet in MDL V2000 (code 4) -> ChiralCraft internal 1
        }
        mol.add_atom(atom);
    }

    for j in 0..n_bonds {
        let line = lines[4 + n_atoms + j];
        let a1: usize = line
            .get(0..3)
            .and_then(|s| s.trim().parse::<usize>().ok())
            .ok_or_else(|| format!("Invalid bond atom1 at bond {}", j))?;
        let a2: usize = line
            .get(3..6)
            .and_then(|s| s.trim().parse::<usize>().ok())
            .ok_or_else(|| format!("Invalid bond atom2 at bond {}", j))?;
        if a1 == 0 || a2 == 0 || a1 > n_atoms || a2 > n_atoms {
            return Err(format!("Bond endpoint out of range at bond {}", j));
        }
        let ord_code: i32 = line.get(6..9).and_then(|s| s.trim().parse().ok()).unwrap_or(1);
        let stereo_code: i32 = line.get(9..12).and_then(|s| s.trim().parse().ok()).unwrap_or(0);
        let order = match ord_code {
            1..=3 => ord_code,
            4 => 6, // Aromatic
            _ => 1,
        };
        let display = match stereo_code {
            1 => BondDisplay::WedgeBegin,
            4 => BondDisplay::Wavy,
            6 => BondDisplay::WedgedHashBegin,
            _ => BondDisplay::Solid,
        };
        mol.add_bond(j as u32 + 1, a1 - 1, a2 - 1, order, display);
    }

    // Parse M  CHG / M  ISO / M  RAD property block lines
    for line in &lines[4 + n_atoms + n_bonds..] {
        if line.starts_with("M  END") {
            break;
        }
        if line.starts_with("M  CHG") || line.starts_with("M  ISO") || line.starts_with("M  RAD") {
            let tokens: Vec<&str> = line[6..].split_whitespace().collect();
            if let Some(cnt_s) = tokens.first() {
                if let Ok(cnt) = cnt_s.parse::<usize>() {
                    for k in 0..cnt {
                        if let (Some(idx_s), Some(val_s)) = (tokens.get(1 + k * 2), tokens.get(2 + k * 2)) {
                            if let (Ok(idx1), Ok(v)) = (idx_s.parse::<usize>(), val_s.parse::<i32>()) {
                                if idx1 >= 1 && idx1 <= mol.atoms.len() {
                                    if line.starts_with("M  CHG") {
                                        mol.atoms[idx1 - 1].charge = v;
                                    } else if line.starts_with("M  ISO") {
                                        mol.atoms[idx1 - 1].isotope = v;
                                    } else if line.starts_with("M  RAD") {
                                        mol.atoms[idx1 - 1].radical = mdl_rad_to_chem(v);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(mol)
}

/// Parsed MDL V2000 Data SGroup (`M  STY`/`M  SAL`/`M  SDT`/`M  SDD`/`M  SED`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileDataSGroup {
    pub sgroup_index: usize,
    pub atom_indices_1based: Vec<usize>,
    pub field_name: String,
    pub tag: Option<String>,
    pub data: String,
}

/// Extracts MDL V2000 Data SGroups (`DAT`) produced by ChiralCraft `Object > Attach Data...` -> `.mol`.
pub fn molfile_v2000_extract_data_sgroups(molfile: &str) -> Vec<MolfileDataSGroup> {
    let mut order: Vec<usize> = Vec::new();
    let mut map: HashMap<usize, MolfileDataSGroup> = HashMap::new();

    for line in molfile.lines() {
        if line.starts_with("M  STY") {
            let toks: Vec<&str> = line[6..].split_whitespace().collect();
            if let Some(cnt) = toks.first().and_then(|s| s.parse::<usize>().ok()) {
                for k in 0..cnt {
                    if let (Some(idx_s), Some("DAT")) = (toks.get(1 + k * 2), toks.get(2 + k * 2).copied()) {
                        if let Ok(idx) = idx_s.parse::<usize>() {
                            if !map.contains_key(&idx) {
                                order.push(idx);
                                map.insert(
                                    idx,
                                    MolfileDataSGroup {
                                        sgroup_index: idx,
                                        atom_indices_1based: Vec::new(),
                                        field_name: String::new(),
                                        tag: None,
                                        data: String::new(),
                                    },
                                );
                            }
                        }
                    }
                }
            }
        } else if line.starts_with("M  SAL") {
            let toks: Vec<&str> = line[6..].split_whitespace().collect();
            if let (Some(Ok(idx)), Some(Ok(cnt))) = (toks.first().map(|s| s.parse::<usize>()), toks.get(1).map(|s| s.parse::<usize>())) {
                if let Some(sg) = map.get_mut(&idx) {
                    for t in toks.iter().skip(2).take(cnt) {
                        if let Ok(a1) = t.parse::<usize>() {
                            sg.atom_indices_1based.push(a1);
                        }
                    }
                }
            }
        } else if line.starts_with("M  SDT") && line.len() >= 11 {
            if let Ok(idx) = line[6..10].trim().parse::<usize>() {
                if let Some(sg) = map.get_mut(&idx) {
                    sg.field_name = line[10..].trim().to_string();
                }
            }
        } else if line.starts_with("M  SDD") && line.len() >= 11 {
            let toks: Vec<&str> = line[6..].split_whitespace().collect();
            if let Some(Ok(idx)) = toks.first().map(|s| s.parse::<usize>()) {
                if let Some(sg) = map.get_mut(&idx) {
                    if toks.len() >= 8 {
                        sg.tag = Some(toks[6].to_string());
                    }
                }
            }
        } else if line.starts_with("M  SED") && line.len() >= 11 {
            if let Ok(idx) = line[6..10].trim().parse::<usize>() {
                if let Some(sg) = map.get_mut(&idx) {
                    sg.data = line[10..].trim().to_string();
                }
            }
        }
    }

    order.into_iter().filter_map(|i| map.remove(&i)).collect()
}

/// Serializes a `ChemMolecule` into an MDL Molfile V2000 string.
pub fn molecule_to_molfile_v2000(mol: &ChemMolecule, title: &str) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str("  ChiralCraft 1.0 2D\n\n");
    out.push_str(&format!("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n", mol.atoms.len(), mol.bonds.len()));

    for atom in &mol.atoms {
        let sym = if let Some(el) = chiralcraft_chem::get_element(atom.z) {
            el.symbol
        } else if atom.z == 225 {
            "D"
        } else if atom.z == 226 {
            "T"
        } else {
            "*"
        };
        let chg_code = match atom.charge {
            3 => 1,
            2 => 2,
            1 => 3,
            -1 => 5,
            -2 => 6,
            -3 => 7,
            _ => 0,
        };
        out.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:<3} 0{:3}  0  0  0  0  0  0  0  0  0  0\n",
            atom.x, -atom.y, atom.z_coord, sym, chg_code
        ));
    }

    for bond in &mol.bonds {
        let ord_code = match bond.order {
            1..=3 => bond.order,
            6 => 4,
            _ => 1,
        };
        let stereo_code = match bond.display {
            BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => 1,
            BondDisplay::Wavy => 4,
            BondDisplay::WedgedHashBegin | BondDisplay::Hash => 6,
            _ => 0,
        };
        out.push_str(&format!(
            "{:3}{:3}{:3}{:3}  0  0  0\n",
            bond.begin_atom + 1,
            bond.end_atom + 1,
            ord_code,
            stereo_code
        ));
    }

    // Emit M  CHG and M  ISO lines if any non-zero
    let charges: Vec<(usize, i32)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.charge != 0)
        .map(|(i, a)| (i + 1, a.charge))
        .collect();
    for chunk in charges.chunks(8) {
        out.push_str(&format!("M  CHG{:3}", chunk.len()));
        for &(idx, chg) in chunk {
            out.push_str(&format!("{:4}{:4}", idx, chg));
        }
        out.push('\n');
    }
    let isotopes: Vec<(usize, i32)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.isotope != 0)
        .map(|(i, a)| (i + 1, a.isotope))
        .collect();
    for chunk in isotopes.chunks(8) {
        out.push_str(&format!("M  ISO{:3}", chunk.len()));
        for &(idx, iso) in chunk {
            out.push_str(&format!("{:4}{:4}", idx, iso));
        }
        out.push('\n');
    }
    let radicals: Vec<(usize, i32)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.radical != 0)
        .map(|(i, a)| (i + 1, chem_rad_to_mdl(a.radical)))
        .filter(|&(_, r)| r != 0)
        .collect();
    for chunk in radicals.chunks(8) {
        out.push_str(&format!("M  RAD{:3}", chunk.len()));
        for &(idx, rad) in chunk {
            out.push_str(&format!("{:4}{:4}", idx, rad));
        }
        out.push('\n');
    }
    out.push_str("M  END\n");
    out
}


/// Extracted enhanced stereochemistry marker (`S03`/`S04`):
/// `(atom_index_in_fragment, group_type /* 2=Abs, 3=Or, 4=And */, group_number)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnhancedStereoMarker {
    pub atom_idx: usize,
    pub group_type: u32,
    pub group_number: u32,
}

impl CdxmlDocument {
    /// Extracts enhanced stereochemistry attributes (`EnhancedStereoType`, `EnhancedStereoGroupNum`)
    /// per fragment occurrence (`S03`/`S04`, matching `CoreChemistryCommon::AdjustStereoCodeForAtom`).
    pub fn extract_enhanced_stereo(&self) -> HashMap<OccurrenceId, Vec<EnhancedStereoMarker>> {
        let mut out: HashMap<OccurrenceId, Vec<EnhancedStereoMarker>> = HashMap::new();
        for (occ, el) in self.elements.iter().enumerate() {
            if el.tag != "fragment" {
                continue;
            }
            let mut atom_idx = 0usize;
            let mut markers = Vec::new();
            for c in &el.children {
                if let CdxmlNodeChild::Element(cid) = c {
                    let node = &self.elements[*cid];
                    if node.tag == "n" {
                        if let Some(st_type) = node.attr("EnhancedStereoType") {
                            let group_type = match st_type {
                                "Abs" => 2,
                                "Or" => 3,
                                "And" => 4,
                                _ => 0,
                            };
                            if group_type != 0 {
                                let group_number = node
                                    .attr("EnhancedStereoGroupNum")
                                    .and_then(|v| v.parse::<u32>().ok())
                                    .unwrap_or(1);
                                markers.push(EnhancedStereoMarker {
                                    atom_idx,
                                    group_type,
                                    group_number,
                                });
                            }
                        }
                        atom_idx += 1;
                    }
                }
            }
            if !markers.is_empty() {
                out.insert(occ, markers);
            }
        }
        out
    }
}

fn collect_v3000_lines(text: &str) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    let mut cont = String::new();
    for raw_line in text.lines() {
        let s = raw_line.strip_prefix("M  V30 ").unwrap_or("");
        if s.is_empty() {
            continue;
        }
        if let Some(body) = s.strip_suffix('-') {
            cont.push_str(body);
        } else {
            cont.push_str(s);
            lines.push(std::mem::take(&mut cont));
        }
    }
    if !cont.is_empty() {
        lines.push(cont);
    }
    lines
}

fn parse_v3000_ctab_lines(lines: &[String]) -> Result<(ChemMolecule, Vec<u32>, Vec<u32>), String> {
    let mut mol = ChemMolecule::new();
    let mut atom_maps: Vec<u32> = Vec::new();
    let mut bond_rxctrs: Vec<u32> = Vec::new();
    let mut in_atom = false;
    let mut in_bond = false;
    let mut id_map: HashMap<u32, usize> = HashMap::new();

    for line in lines {
        let trim = line.trim();
        if trim == "BEGIN ATOM" {
            in_atom = true;
            continue;
        }
        if trim == "END ATOM" {
            in_atom = false;
            continue;
        }
        if trim == "BEGIN BOND" {
            in_bond = true;
            continue;
        }
        if trim == "END BOND" {
            in_bond = false;
            continue;
        }
        if in_atom {
            let parts: Vec<&str> = trim.split_whitespace().collect();
            if parts.len() < 5 {
                return Err(format!("Truncated V3000 atom line: {trim}"));
            }
            let aid: u32 = parts[0].parse().map_err(|_| format!("Bad V3000 atom id: {}", parts[0]))?;
            let sym = parts[1];
            let z = if sym == "D" || sym == "T" {
                1
            } else {
                get_element_by_symbol(sym).map(|e| e.z).unwrap_or(0)
            };
            let mut atom = ChemAtom::new(aid, z);
            atom.x = parts[2].parse().unwrap_or(0.0);
            atom.y = parts[3].parse().unwrap_or(0.0);
            atom.z_coord = parts[4].parse().unwrap_or(0.0);
            if sym == "D" {
                atom.isotope = 2;
            } else if sym == "T" {
                atom.isotope = 3;
            }
            let mut aam = 0u32;
            let mut extra_start = 5usize;
            if let Some(p5) = parts.get(5) {
                if !p5.contains('=') {
                    aam = p5.parse().unwrap_or(0);
                    extra_start = 6;
                }
            }
            for kv in parts.iter().skip(extra_start) {
                if let Some((k, v)) = kv.split_once('=') {
                    match k {
                        "CHG" => atom.charge = v.parse().unwrap_or(0),
                        "RAD" => atom.radical = mdl_rad_to_chem(v.parse().unwrap_or(0)),
                        "MASS" => atom.isotope = v.parse().unwrap_or(0),
                        "VAL" => {
                            let val = v.parse::<i32>().unwrap_or(0);
                            if val > 0 {
                                atom.valence_override = Some(val);
                            }
                        }
                        _ => {}
                    }
                }
            }
            let idx = mol.add_atom(atom);
            atom_maps.push(aam);
            id_map.insert(aid, idx);
        } else if in_bond {
            let parts: Vec<&str> = trim.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(format!("Truncated V3000 bond line: {trim}"));
            }
            let bid: u32 = parts[0].parse().unwrap_or(1);
            let order_code: i32 = parts[1].parse().unwrap_or(1);
            let order = match order_code {
                4 => 6,
                _ => order_code,
            };
            let a1_id: u32 = parts[2].parse().map_err(|_| format!("Bad V3000 bond begin: {}", parts[2]))?;
            let a2_id: u32 = parts[3].parse().map_err(|_| format!("Bad V3000 bond end: {}", parts[3]))?;
            let a1 = *id_map.get(&a1_id).ok_or_else(|| format!("Missing V3000 atom {a1_id}"))?;
            let a2 = *id_map.get(&a2_id).ok_or_else(|| format!("Missing V3000 atom {a2_id}"))?;
            let mut disp = BondDisplay::Solid;
            let mut rxctr = 0u32;
            for kv in parts.iter().skip(4) {
                if let Some((k, v)) = kv.split_once('=') {
                    match k {
                        "CFG" => {
                            disp = match v {
                                "1" => BondDisplay::WedgeBegin,
                                "2" => BondDisplay::Wavy,
                                "3" => BondDisplay::WedgedHashBegin,
                                _ => BondDisplay::Solid,
                            };
                        }
                        "RXCTR" => rxctr = v.parse().unwrap_or(0),
                        _ => {}
                    }
                }
            }
            mol.add_bond(bid, a1, a2, order, disp);
            bond_rxctrs.push(rxctr);
        }
    }
    Ok((mol, atom_maps, bond_rxctrs))
}

fn write_v3000_ctab_body(
    out: &mut String,
    mol: &ChemMolecule,
    atom_maps: Option<&[u32]>,
    rxn_ctrs: Option<&[u32]>,
) {
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms.iter().enumerate() {
        let sym = get_element(atom.z).map(|e| e.symbol).unwrap_or("*");
        let aam = atom_maps.and_then(|m| m.get(idx).copied()).unwrap_or(0);
        out.push_str(&format!(
            "M  V30 {} {} {:.6} {:.6} {:.6} {}",
            idx + 1,
            sym,
            atom.x,
            atom.y,
            atom.z_coord,
            aam
        ));
        if atom.charge != 0 {
            out.push_str(&format!(" CHG={}", atom.charge));
        }
        if atom.radical != 0 {
            let mdl_r = chem_rad_to_mdl(atom.radical);
            if mdl_r != 0 {
                out.push_str(&format!(" RAD={mdl_r}"));
            }
        }
        if atom.isotope != 0 {
            out.push_str(&format!(" MASS={}", atom.isotope));
        }
        if let Some(vo) = atom.valence_override {
            out.push_str(&format!(" VAL={vo}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (idx, bond) in mol.bonds.iter().enumerate() {
        let ord = match bond.order {
            6 => 4,
            1..=3 => bond.order,
            _ => 1,
        };
        out.push_str(&format!(
            "M  V30 {} {} {} {}",
            idx + 1,
            ord,
            bond.begin_atom + 1,
            bond.end_atom + 1
        ));
        match bond.display {
            BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => out.push_str(" CFG=1"),
            BondDisplay::Wavy => out.push_str(" CFG=2"),
            BondDisplay::WedgedHashBegin | BondDisplay::Hash => out.push_str(" CFG=3"),
            _ => {}
        }
        if let Some(rc) = rxn_ctrs.and_then(|r| r.get(idx).copied()) {
            if rc != 0 {
                out.push_str(&format!(" RXCTR={rc}"));
            }
        }
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    out.push_str("M  V30 END CTAB\n");
}

/// Parses an MDL Molfile V3000 (`M  V30 ...`) into a `ChemMolecule`.
pub fn molfile_v3000_to_molecule(molfile: &str) -> Result<ChemMolecule, String> {
    let lines = collect_v3000_lines(molfile);
    let (mol, _, _) = parse_v3000_ctab_lines(&lines)?;
    Ok(mol)
}

/// Serializes a `ChemMolecule` into MDL Molfile V3000 (`CTAB`) format.
pub fn molecule_to_molfile_v3000(mol: &ChemMolecule, title: &str) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str("  ChiralCraft 2D         \n\n");
    out.push_str("  0  0  0     0  0            999 V3000\n");
    write_v3000_ctab_body(&mut out, mol, None, None);
    out.push_str("M  END\n");
    out
}


/// Auto-detects MDL Molfile V2000 vs V3000 (`M  V30 BEGIN CTAB`).
pub fn molfile_to_molecule(molfile: &str) -> Result<ChemMolecule, String> {
    if molfile.contains("V3000") || molfile.contains("M  V30") {
        molfile_v3000_to_molecule(molfile)
    } else {
        molfile_v2000_to_molecule(molfile)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub title: String,
    pub molecule: ChemMolecule,
    pub properties: BTreeMap<String, String>,
}

/// Parses an MDL SDfile (`$$$$` delimited records with `> <KEY>` data tags).
pub fn sdf_to_molecules(sdf: &str) -> Result<Vec<SdfRecord>, String> {
    let mut records = Vec::new();
    for raw_block in sdf.split("$$$$") {
        let block = raw_block.trim_start_matches(['\r', '\n']);
        if block.trim().is_empty() {
            continue;
        }
        let title = block.lines().next().unwrap_or("").trim().to_string();
        let (ctab_part, data_part) = if let Some(idx) = block.find("M  END") {
            let end_idx = idx + "M  END".len();
            (&block[..end_idx], &block[end_idx..])
        } else {
            (block, "")
        };
        let molecule = molfile_to_molecule(ctab_part)?;
        let mut properties = BTreeMap::new();
        let mut cur_key: Option<String> = None;
        let mut cur_vals: Vec<String> = Vec::new();
        for line in data_part.lines() {
            if let Some(rest) = line.trim().strip_prefix('>') {
                if let Some(prev_key) = cur_key.take() {
                    properties.insert(prev_key, cur_vals.join("\n").trim_end().to_string());
                    cur_vals.clear();
                }
                if let Some(l) = rest.find('<') {
                    if let Some(r) = rest[l + 1..].find('>') {
                        cur_key = Some(rest[l + 1..l + 1 + r].to_string());
                    }
                }
            } else if cur_key.is_some() {
                if line.trim().is_empty() && cur_vals.is_empty() {
                    continue;
                }
                if line.trim().is_empty() {
                    if let Some(prev_key) = cur_key.take() {
                        properties.insert(prev_key, cur_vals.join("\n").trim_end().to_string());
                        cur_vals.clear();
                    }
                } else {
                    cur_vals.push(line.to_string());
                }
            }
        }
        if let Some(prev_key) = cur_key {
            properties.insert(prev_key, cur_vals.join("\n").trim_end().to_string());
        }
        records.push(SdfRecord {
            title,
            molecule,
            properties,
        });
    }
    Ok(records)
}

/// Serializes multiple `SdfRecord` entries into an MDL SDfile (`$$$$` delimited).
pub fn molecules_to_sdf(records: &[SdfRecord]) -> String {
    let mut out = String::new();
    for rec in records {
        out.push_str(&molecule_to_molfile_v2000(&rec.molecule, &rec.title));
        for (k, v) in &rec.properties {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
        out.push_str("$$$$\n");
    }
    out
}


#[derive(Debug, Clone, PartialEq)]
pub struct MdlRxnFile {
    pub title: String,
    pub reactants: Vec<ChemMolecule>,
    pub products: Vec<ChemMolecule>,
}

/// Parses an MDL RXN V2000 file (`$RXN` header + `$MOL` blocks) into reactants and products (`MDLRxnFile`).
pub fn rxn_v2000_to_reaction(rxn: &str) -> Result<MdlRxnFile, String> {
    let lines: Vec<&str> = rxn.lines().collect();
    if lines.len() < 5 || !lines[0].trim_start().starts_with("$RXN") {
        return Err("Missing $RXN header".to_string());
    }
    let title = lines[1].trim().to_string();
    let counts = lines[4];
    let n_react: usize = counts.get(0..3).and_then(|s| s.trim().parse().ok()).unwrap_or(0);
    let n_prod: usize = counts.get(3..6).and_then(|s| s.trim().parse().ok()).unwrap_or(0);

    let rest = lines[5..].join("\n");
    let mut mol_blocks: Vec<ChemMolecule> = Vec::new();
    for part in rest.split("$MOL").skip(1) {
        let block = part
            .strip_prefix("\r\n")
            .or_else(|| part.strip_prefix('\n'))
            .unwrap_or(part);
        mol_blocks.push(molfile_to_molecule(block)?);
    }
    if mol_blocks.len() < n_react + n_prod {
        return Err(format!(
            "RXN declared {} reactants + {} products, found {} $MOL blocks",
            n_react,
            n_prod,
            mol_blocks.len()
        ));
    }
    let reactants = mol_blocks[..n_react].to_vec();
    let products = mol_blocks[n_react..n_react + n_prod].to_vec();
    Ok(MdlRxnFile {
        title,
        reactants,
        products,
    })
}

/// Serializes an `MdlRxnFile` (`reactants` + `products`) to MDL RXN V2000 (`$RXN` + `$MOL`).
pub fn reaction_to_rxn_v2000(rxn: &MdlRxnFile) -> String {
    let mut out = String::new();
    out.push_str("$RXN\n");
    out.push_str(&format!("{}\n", rxn.title));
    out.push_str("  ChiralCraft\n\n");
    out.push_str(&format!("{:3}{:3}\n", rxn.reactants.len(), rxn.products.len()));
    for r in &rxn.reactants {
        out.push_str("$MOL\n");
        out.push_str(&molecule_to_molfile_v2000(r, ""));
    }
    for p in &rxn.products {
        out.push_str("$MOL\n");
        out.push_str(&molecule_to_molfile_v2000(p, ""));
    }
    out
}

/// Extended MDL RXN V3000 representation preserving atom-mapping numbers (`aam`) and `RXCTR` flags.
#[derive(Debug, Clone, PartialEq)]
pub struct MdlRxnV3000File {
    pub base: MdlRxnFile,
    pub reactant_atom_maps: Vec<Vec<u32>>,
    pub product_atom_maps: Vec<Vec<u32>>,
    pub reactant_bond_rxctrs: Vec<Vec<u32>>,
    pub product_bond_rxctrs: Vec<Vec<u32>>,
}

/// Parses an MDL RXN V3000 (`$RXN V3000`) string into reactants, products, atom maps, and `RXCTR` values.
pub fn rxn_v3000_to_reaction(rxn: &str) -> Result<MdlRxnV3000File, String> {
    let raw_lines: Vec<&str> = rxn.lines().collect();
    if raw_lines.is_empty() || !raw_lines[0].trim_start().starts_with("$RXN V3000") {
        return Err("Missing $RXN V3000 header".to_string());
    }
    let title = raw_lines.get(1).map(|s| s.trim()).unwrap_or("").to_string();
    let v30_lines = collect_v3000_lines(rxn);
    let mut section = 0u8; // 1=reactant, 2=product
    let mut in_ctab = false;
    let mut ctab_buf: Vec<String> = Vec::new();
    let mut reactants = Vec::new();
    let mut products = Vec::new();
    let mut reactant_atom_maps = Vec::new();
    let mut product_atom_maps = Vec::new();
    let mut reactant_bond_rxctrs = Vec::new();
    let mut product_bond_rxctrs = Vec::new();

    for line in v30_lines {
        let trim = line.trim();
        if trim == "BEGIN REACTANT" {
            section = 1;
            continue;
        } else if trim == "END REACTANT" {
            section = 0;
            continue;
        } else if trim == "BEGIN PRODUCT" {
            section = 2;
            continue;
        } else if trim == "END PRODUCT" {
            section = 0;
            continue;
        }
        if trim == "BEGIN CTAB" {
            in_ctab = true;
            ctab_buf.clear();
            continue;
        }
        if trim == "END CTAB" && in_ctab {
            in_ctab = false;
            let (mol, aam, rxctrs) = parse_v3000_ctab_lines(&ctab_buf)?;
            if section == 1 {
                reactants.push(mol);
                reactant_atom_maps.push(aam);
                reactant_bond_rxctrs.push(rxctrs);
            } else if section == 2 {
                products.push(mol);
                product_atom_maps.push(aam);
                product_bond_rxctrs.push(rxctrs);
            }
            continue;
        }
        if in_ctab {
            ctab_buf.push(line);
        }
    }

    Ok(MdlRxnV3000File {
        base: MdlRxnFile {
            title,
            reactants,
            products,
        },
        reactant_atom_maps,
        product_atom_maps,
        reactant_bond_rxctrs,
        product_bond_rxctrs,
    })
}

/// Serializes an `MdlRxnV3000File` (`$RXN V3000`) with atom mapping numbers and bond `RXCTR` flags.
pub fn reaction_to_rxn_v3000(rxn: &MdlRxnV3000File) -> String {
    let mut out = String::new();
    out.push_str("$RXN V3000\n");
    out.push_str(&format!("{}\n", rxn.base.title));
    out.push_str("  ChiralCraft 2D\n\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {}\n",
        rxn.base.reactants.len(),
        rxn.base.products.len()
    ));
    out.push_str("M  V30 BEGIN REACTANT\n");
    for (i, r) in rxn.base.reactants.iter().enumerate() {
        write_v3000_ctab_body(
            &mut out,
            r,
            rxn.reactant_atom_maps.get(i).map(|v| v.as_slice()),
            rxn.reactant_bond_rxctrs.get(i).map(|v| v.as_slice()),
        );
    }
    out.push_str("M  V30 END REACTANT\n");
    out.push_str("M  V30 BEGIN PRODUCT\n");
    for (i, p) in rxn.base.products.iter().enumerate() {
        write_v3000_ctab_body(
            &mut out,
            p,
            rxn.product_atom_maps.get(i).map(|v| v.as_slice()),
            rxn.product_bond_rxctrs.get(i).map(|v| v.as_slice()),
        );
    }
    out.push_str("M  V30 END PRODUCT\n");
    out.push_str("M  END\n");
    out
}

/// Auto-detects MDL RXN V2000 vs `$RXN V3000` into `MdlRxnFile`.
pub fn rxn_to_reaction(rxn: &str) -> Result<MdlRxnFile, String> {
    if rxn.trim_start().starts_with("$RXN V3000") {
        Ok(rxn_v3000_to_reaction(rxn)?.base)
    } else {
        rxn_v2000_to_reaction(rxn)
    }
}

/// Parses a SMILES string into a `ChemMolecule` (supports organic/aromatic subsets, brackets with
/// isotopes/charges/explicit H, single/double/triple/aromatic bonds, branches, rings `1..=99`, `.` fragments).
pub fn smiles_to_molecule(smiles: &str) -> Result<ChemMolecule, String> {
    let trimmed = smiles.trim();
    if trimmed.is_empty() {
        return Err("Empty SMILES string".to_string());
    }
    let bytes = trimmed.as_bytes();
    let mut i = 0usize;
    let mut mol = ChemMolecule::new();
    let mut next_atom_id = 1u32;
    let mut next_bond_id = 1u32;
    let mut prev_atom: Option<usize> = None;
    let mut branch_stack: Vec<Option<usize>> = Vec::new();
    let mut ring_opens: HashMap<usize, (usize, Option<i32>)> = HashMap::new();
    let mut pending_bond: Option<i32> = None;
    let mut aromatic_atoms: Vec<bool> = Vec::new();

    while i < bytes.len() {
        match bytes[i] {
            b'.' => {
                prev_atom = None;
                pending_bond = None;
                i += 1;
            }
            b'(' => {
                branch_stack.push(prev_atom);
                pending_bond = None;
                i += 1;
            }
            b')' => {
                prev_atom = branch_stack
                    .pop()
                    .ok_or_else(|| "Unmatched ')' in SMILES".to_string())?;
                pending_bond = None;
                i += 1;
            }
            b'-' => {
                pending_bond = Some(1);
                i += 1;
            }
            b'=' => {
                pending_bond = Some(2);
                i += 1;
            }
            b'#' => {
                pending_bond = Some(3);
                i += 1;
            }
            b':' => {
                pending_bond = Some(6);
                i += 1;
            }
            b'/' | b'\\' => {
                if pending_bond.is_none() {
                    pending_bond = Some(1);
                }
                i += 1;
            }
            b'%' | b'0'..=b'9' => {
                let rnum = if bytes[i] == b'%' {
                    i += 1;
                    if i + 2 > bytes.len() || !bytes[i].is_ascii_digit() || !bytes[i + 1].is_ascii_digit() {
                        return Err("Invalid %nn ring closure in SMILES".to_string());
                    }
                    let n = ((bytes[i] - b'0') as usize) * 10 + ((bytes[i + 1] - b'0') as usize);
                    i += 2;
                    n
                } else {
                    let n = (bytes[i] - b'0') as usize;
                    i += 1;
                    n
                };
                let cur = prev_atom.ok_or_else(|| "Ring digit without preceding atom".to_string())?;
                if let Some((open_atom, open_order)) = ring_opens.remove(&rnum) {
                    let order = pending_bond.or(open_order).unwrap_or_else(|| {
                        if aromatic_atoms[open_atom] && aromatic_atoms[cur] {
                            6
                        } else {
                            1
                        }
                    });
                    mol.add_bond(next_bond_id, open_atom, cur, order, BondDisplay::Solid);
                    next_bond_id += 1;
                    pending_bond = None;
                } else {
                    ring_opens.insert(rnum, (cur, pending_bond));
                    pending_bond = None;
                }
            }
            b'[' => {
                i += 1;
                let close = trimmed[i..]
                    .find(']')
                    .ok_or_else(|| "Unterminated '[' atom in SMILES".to_string())?;
                let inside = &trimmed[i..i + close];
                i += close + 1;
                let mut p = 0usize;
                let ib = inside.as_bytes();
                let mut iso = 0i32;
                while p < ib.len() && ib[p].is_ascii_digit() {
                    iso = iso * 10 + (ib[p] - b'0') as i32;
                    p += 1;
                }
                let sym_start = p;
                if p >= ib.len() {
                    return Err("Empty SMILES bracket atom".to_string());
                }
                let is_arom = ib[p].is_ascii_lowercase();
                p += 1;
                if p < ib.len()
                    && ib[p].is_ascii_lowercase()
                    && !matches!(ib[p], b'h')
                {
                    let candidate = &inside[sym_start..=p];
                    let cap = format!(
                        "{}{}",
                        candidate.chars().next().unwrap().to_ascii_uppercase(),
                        &candidate[1..]
                    );
                    if get_element_by_symbol(&cap).is_some() {
                        p += 1;
                    }
                }
                let raw_sym = &inside[sym_start..p];
                let lookup = if is_arom {
                    raw_sym.to_ascii_uppercase()
                } else {
                    raw_sym.to_string()
                };
                let z = get_element_by_symbol(&lookup)
                    .map(|e| e.z)
                    .ok_or_else(|| format!("Unknown element in SMILES bracket: {inside}"))?;
                let mut chg = 0i32;
                while p < ib.len() {
                    match ib[p] {
                        b'+' => {
                            p += 1;
                            if p < ib.len() && ib[p].is_ascii_digit() {
                                let mut mag = 0i32;
                                while p < ib.len() && ib[p].is_ascii_digit() {
                                    mag = mag * 10 + (ib[p] - b'0') as i32;
                                    p += 1;
                                }
                                chg += mag;
                            } else {
                                chg += 1;
                            }
                        }
                        b'-' => {
                            p += 1;
                            if p < ib.len() && ib[p].is_ascii_digit() {
                                let mut mag = 0i32;
                                while p < ib.len() && ib[p].is_ascii_digit() {
                                    mag = mag * 10 + (ib[p] - b'0') as i32;
                                    p += 1;
                                }
                                chg -= mag;
                            } else {
                                chg -= 1;
                            }
                        }
                        _ => p += 1,
                    }
                }
                let mut atom = ChemAtom::new(next_atom_id, z);
                next_atom_id += 1;
                atom.isotope = iso;
                atom.charge = chg;
                let cur_idx = mol.add_atom(atom);
                aromatic_atoms.push(is_arom);
                if let Some(prev) = prev_atom {
                    let order = pending_bond.unwrap_or_else(|| {
                        if aromatic_atoms[prev] && is_arom { 6 } else { 1 }
                    });
                    mol.add_bond(next_bond_id, prev, cur_idx, order, BondDisplay::Solid);
                    next_bond_id += 1;
                }
                pending_bond = None;
                prev_atom = Some(cur_idx);
            }
            ch => {
                let (sym, is_arom, adv) = match ch {
                    b'B' if i + 1 < bytes.len() && bytes[i + 1] == b'r' => ("Br", false, 2),
                    b'C' if i + 1 < bytes.len() && bytes[i + 1] == b'l' => ("Cl", false, 2),
                    b'B' => ("B", false, 1),
                    b'C' => ("C", false, 1),
                    b'N' => ("N", false, 1),
                    b'O' => ("O", false, 1),
                    b'P' => ("P", false, 1),
                    b'S' => ("S", false, 1),
                    b'F' => ("F", false, 1),
                    b'I' => ("I", false, 1),
                    b'b' => ("B", true, 1),
                    b'c' => ("C", true, 1),
                    b'n' => ("N", true, 1),
                    b'o' => ("O", true, 1),
                    b'p' => ("P", true, 1),
                    b's' => ("S", true, 1),
                    _ => return Err(format!("Unsupported SMILES token '{}' at index {i}", ch as char)),
                };
                i += adv;
                let z = get_element_by_symbol(sym).unwrap().z;
                let atom = ChemAtom::new(next_atom_id, z);
                next_atom_id += 1;
                let cur_idx = mol.add_atom(atom);
                aromatic_atoms.push(is_arom);
                if let Some(prev) = prev_atom {
                    let order = pending_bond.unwrap_or_else(|| {
                        if aromatic_atoms[prev] && is_arom { 6 } else { 1 }
                    });
                    mol.add_bond(next_bond_id, prev, cur_idx, order, BondDisplay::Solid);
                    next_bond_id += 1;
                }
                pending_bond = None;
                prev_atom = Some(cur_idx);
            }
        }
    }

    if !branch_stack.is_empty() {
        return Err("Unclosed '(' branch in SMILES".to_string());
    }
    if !ring_opens.is_empty() {
        return Err("Unclosed ring closure digit in SMILES".to_string());
    }
    Ok(mol)
}

/// Serializes a `ChemMolecule` to SMILES using `CcAromaticityInfo` for aromatic rings and
/// DFS spanning tree with numbered ring closures.
pub fn molecule_to_smiles(mol: &ChemMolecule) -> String {
    if mol.atoms.is_empty() {
        return String::new();
    }
    let arom = CcAromaticityInfo::perceive(mol);
    let val = CcVal::compute(mol);
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    let atom_is_arom: Vec<bool> = (0..n)
        .map(|i| adj[i].iter().any(|(_, bi)| arom.aromatic_bonds.contains(bi)))
        .collect();

    let format_atom = |i: usize| -> String {
        let a = &mol.atoms[i];
        let sym = get_element(a.z).map(|e| e.symbol).unwrap_or("*");
        let is_organic = matches!(a.z, 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53);
        let is_arom = atom_is_arom[i];
        let disp_sym = if is_arom && matches!(a.z, 5 | 6 | 7 | 8 | 15 | 16) {
            sym.to_ascii_lowercase()
        } else {
            sym.to_string()
        };
        let needs_h_in_arom = is_arom && a.z == 7 && val.implicit_h_counts[i] > 0;
        if is_organic && a.charge == 0 && a.isotope == 0 && a.radical == 0 && !needs_h_in_arom {
            disp_sym
        } else {
            let mut s = String::from("[");
            if a.isotope != 0 {
                s.push_str(&a.isotope.to_string());
            }
            s.push_str(&disp_sym);
            if val.implicit_h_counts[i] > 0 {
                s.push('H');
                if val.implicit_h_counts[i] > 1 {
                    s.push_str(&val.implicit_h_counts[i].to_string());
                }
            }
            if a.charge > 0 {
                s.push('+');
                if a.charge > 1 {
                    s.push_str(&a.charge.to_string());
                }
            } else if a.charge < 0 {
                s.push('-');
                if a.charge < -1 {
                    s.push_str(&(-a.charge).to_string());
                }
            }
            s.push(']');
            s
        }
    };

    let format_bond = |bi: usize| -> &'static str {
        if arom.aromatic_bonds.contains(&bi) {
            return "";
        }
        match mol.bonds[bi].order {
            2 => "=",
            3 => "#",
            _ => "",
        }
    };

    // Pass 1: find tree vs back-edges so ring closures get deterministic digits 1..N
    let mut visited = vec![false; n];
    let mut tree_adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    let mut back_edges_at: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); n];
    let mut next_ring = 1usize;
    let mut roots = Vec::new();

    fn build_tree(
        u: usize,
        parent: Option<usize>,
        adj: &[Vec<(usize, usize)>],
        visited: &mut [bool],
        seen_edge: &mut [bool],
        tree_adj: &mut [Vec<(usize, usize)>],
        back_edges_at: &mut [Vec<(usize, usize, usize)>],
        next_ring: &mut usize,
    ) {
        visited[u] = true;
        for &(v, bi) in &adj[u] {
            if Some(v) == parent || seen_edge[bi] {
                continue;
            }
            seen_edge[bi] = true;
            if !visited[v] {
                tree_adj[u].push((v, bi));
                build_tree(v, Some(u), adj, visited, seen_edge, tree_adj, back_edges_at, next_ring);
            } else {
                let r = *next_ring;
                *next_ring += 1;
                back_edges_at[u].push((r, bi, v));
                back_edges_at[v].push((r, bi, u));
            }
        }
    }

    let mut seen_edge = vec![false; mol.bonds.len()];
    for i in 0..n {
        if !visited[i] {
            roots.push(i);
            build_tree(
                i,
                None,
                &adj,
                &mut visited,
                &mut seen_edge,
                &mut tree_adj,
                &mut back_edges_at,
                &mut next_ring,
            );
        }
    }

    let mut emitted_ring: HashMap<usize, bool> = HashMap::new();

    fn write_node(
        u: usize,
        format_atom: &dyn Fn(usize) -> String,
        format_bond: &dyn Fn(usize) -> &'static str,
        tree_adj: &[Vec<(usize, usize)>],
        back_edges_at: &[Vec<(usize, usize, usize)>],
        emitted_ring: &mut HashMap<usize, bool>,
        out: &mut String,
    ) {
        out.push_str(&format_atom(u));
        for &(rnum, bi, _) in &back_edges_at[u] {
            if !emitted_ring.contains_key(&rnum) {
                emitted_ring.insert(rnum, true);
                out.push_str(format_bond(bi));
            }
            if rnum < 10 {
                out.push_str(&rnum.to_string());
            } else {
                out.push_str(&format!("%{rnum}"));
            }
        }
        let children = &tree_adj[u];
        for (idx, &(v, bi)) in children.iter().enumerate() {
            let is_last = idx + 1 == children.len();
            if !is_last {
                out.push('(');
            }
            out.push_str(format_bond(bi));
            write_node(v, format_atom, format_bond, tree_adj, back_edges_at, emitted_ring, out);
            if !is_last {
                out.push(')');
            }
        }
    }

    let mut parts = Vec::new();
    for root in roots {
        let mut part = String::new();
        write_node(
            root,
            &format_atom,
            &format_bond,
            &tree_adj,
            &back_edges_at,
            &mut emitted_ring,
            &mut part,
        );
        parts.push(part);
    }
    parts.join(".")
}

/// Helper to check if a single-token caption is chemically interpreted by ChiralCraft (`InterpretChemically="yes"`).
/// Observed in P01 `seed-v2.cdxml`:
/// - `"A"` -> `NodeType="GenericNickname" GenericNickname="A"`
/// - `"B"` -> `Element="5" Warning="An atom in this label has an invalid valence."`
pub fn classify_single_caption_token(text: &str) -> Option<(Option<i32>, bool)> {
    if let Some(el) = get_element_by_symbol(text) {
        Some((Some(el.z), true))
    } else if matches!(text, "A" | "M" | "Q" | "X" | "R") {
        Some((None, false))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_cdx_synthetic_and_malformed_guards() {
        let mut header28 = Vec::from(&CDX_MAGIC_PREFIX[..]);
        header28.extend_from_slice(&[0u8; 16]);
        let mut data = header28.clone();
        data.extend_from_slice(&0x8000u16.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());
        // Extended-length property 0x4567 of 65536 bytes
        data.extend_from_slice(&0x4567u16.to_le_bytes());
        data.extend_from_slice(&0xffffu16.to_le_bytes());
        data.extend_from_slice(&65536u32.to_le_bytes());
        data.extend_from_slice(&vec![b'x'; 65536]);
        data.extend_from_slice(&0u16.to_le_bytes()); // end object
        data.extend_from_slice(&0u16.to_le_bytes()); // end file

        let doc = CdxDocument::parse(&data).unwrap();
        assert_eq!(doc.header_bytes.len(), 28);
        assert_eq!(doc.total_objects, 1);
        assert_eq!(doc.total_properties, 1);
        assert_eq!(doc.to_bytes(), data);

        // Reject truncated, trailing, and huge length
        assert!(CdxDocument::parse(&data[..data.len() - 1]).is_err());
        assert!(CdxDocument::parse(&data[..data.len() - 4]).is_err());
        let mut trailing = data.clone();
        trailing.push(b'x');
        assert!(CdxDocument::parse(&trailing).is_err());
    }

    #[test]
    fn test_all_50_shipped_cdx_files_byte_exact_roundtrip() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let inv_path = root.join("research/formats/inventory.json");
        if !inv_path.exists() {
            return;
        }
        let inv_text = fs::read_to_string(inv_path).unwrap();
        let mut cdx_paths = Vec::new();
        for line in inv_text.lines() {
            if line.contains("\"kind\": \"CDX\"") {
                // previous or next lines have path
            }
        }
        // Simple JSON extraction of (path, kind) pairs from inventory.json
        let mut current_path = String::new();
        for line in inv_text.lines() {
            let s = line.trim();
            if let Some(rest) = s.strip_prefix("\"path\": \"") {
                if let Some(end) = rest.find('"') {
                    current_path = rest[..end].to_string();
                }
            } else if s.starts_with("\"kind\": \"CDX\"") {
                cdx_paths.push(current_path.clone());
            }
        }
        assert_eq!(cdx_paths.len(), 50);

        let mut total_objects = 0;
        let mut total_props = 0;
        let mut total_obj_terms = 0;
        let mut total_cdx_fragments = 0;
        for rel in &cdx_paths {
            let bytes = fs::read(root.join(rel)).unwrap();
            let doc = CdxDocument::parse(&bytes).unwrap_or_else(|e| panic!("Failed {}: {:?}", rel, e));
            assert_eq!(doc.header_bytes.len(), 22);
            total_objects += doc.total_objects;
            total_props += doc.total_properties;
            total_obj_terms += doc.total_object_terminators;
            let angles = doc.chain_angles_degrees();
            assert_eq!(angles, vec![120.0]);
            let frags = doc.extract_molecules();
            total_cdx_fragments += frags.len();
            // Verify CDX -> CDXML conversion produces valid parsable CDXML with matching fragment count
            let cdxml_doc = doc.to_cdxml();
            let cdxml_str = cdxml_doc.to_xml_string();
            let reparsed = CdxmlDocument::parse(&cdxml_str).unwrap();
            assert_eq!(reparsed.extract_molecules().len(), frags.len());
            // Verify 1:1 byte-identical serialization roundtrip
            let roundtripped = doc.to_bytes();
            assert_eq!(roundtripped, bytes, "Roundtrip mismatch on {}", rel);
        }
        assert_eq!(total_objects, 60_711);
        assert_eq!(total_props, 402_315);
        assert_eq!(total_obj_terms, 60_711);
        assert_eq!(total_cdx_fragments, 480);
    }

    #[test]
    fn test_all_655_shipped_cdxml_files_semantic_and_roundtrip() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let inv_path = root.join("research/formats/inventory.json");
        if !inv_path.exists() {
            return;
        }
        let inv_text = fs::read_to_string(inv_path).unwrap();
        let mut cdxml_paths = Vec::new();
        let mut current_path = String::new();
        for line in inv_text.lines() {
            let s = line.trim();
            if let Some(rest) = s.strip_prefix("\"path\": \"") {
                if let Some(end) = rest.find('"') {
                    current_path = rest[..end].to_string();
                }
            } else if s.starts_with("\"kind\": \"CDXML\"") {
                cdxml_paths.push(current_path.clone());
            }
        }
        assert_eq!(cdxml_paths.len(), 655);

        let mut docs_with_dup_ids = 0;
        let mut total_dup_groups = 0;
        let mut docs_with_font_collisions = 0;
        let mut total_bond_endpoints = 0;
        let mut total_resolved_endpoints = 0;
        let mut total_schemes = 0;
        let mut total_steps = 0;

        for rel in &cdxml_paths {
            let raw = fs::read_to_string(root.join(rel)).unwrap();
            let doc = CdxmlDocument::parse(&raw).unwrap_or_else(|e| panic!("Failed {}: {}", rel, e));
            let rep = doc.analyze_references();
            if !rep.duplicate_object_id_groups.is_empty() {
                docs_with_dup_ids += 1;
                total_dup_groups += rep.duplicate_object_id_groups.len();
            }
            if !rep.font_object_id_collisions.is_empty() {
                docs_with_font_collisions += 1;
            }
            total_bond_endpoints += rep.bond_endpoints_total;
            total_resolved_endpoints += rep.bond_endpoints_typed_local_resolved;
            assert_eq!(rep.bond_endpoints_typed_local_ambiguous, 0, "Ambiguous endpoint in {}", rel);
            assert_eq!(rep.bond_endpoints_typed_local_dangling, 0, "Dangling endpoint in {}", rel);
            assert_eq!(rep.non_bond_references_unresolved, 0, "Unresolved non-bond ref in {}", rel);

            for rxn in doc.extract_reactions() {
                total_schemes += 1;
                total_steps += rxn.steps.len();
            }

            // Roundtrip XML and verify element count + reference report invariance
            let serialized = doc.to_xml_string();
            let doc2 = CdxmlDocument::parse(&serialized).unwrap();
            assert_eq!(doc.elements.len(), doc2.elements.len(), "Element count changed in {}", rel);
            let rep2 = doc2.analyze_references();
            assert_eq!(rep.bond_endpoints_typed_local_resolved, rep2.bond_endpoints_typed_local_resolved);
            assert_eq!(rep.duplicate_object_id_groups.len(), rep2.duplicate_object_id_groups.len());
        }

        assert_eq!(docs_with_dup_ids, 8);
        assert_eq!(total_dup_groups, 34);
        assert_eq!(docs_with_font_collisions, 74);
        assert_eq!(total_bond_endpoints, 27_234);
        assert_eq!(total_resolved_endpoints, 27_234);
        assert_eq!(total_schemes, 94);
        assert_eq!(total_steps, 371);
    }

    #[test]
    fn test_molfile_v2000_roundtrip_and_chemistry() {
        let mut mol = ChemMolecule::new();
        let mut c1 = ChemAtom::new(1, 6);
        c1.x = 0.0;
        c1.y = 0.0;
        let mut c2 = ChemAtom::new(2, 6);
        c2.x = 14.4;
        c2.y = 0.0;
        let mut o = ChemAtom::new(3, 8);
        o.x = 21.6;
        o.y = -12.47;
        mol.add_atom(c1);
        mol.add_atom(c2);
        mol.add_atom(o);
        mol.add_bond(1, 0, 1, 1, BondDisplay::Solid);
        mol.add_bond(2, 1, 2, 1, BondDisplay::Solid);

        let v2000 = molecule_to_molfile_v2000(&mol, "Ethanol");
        let parsed = molfile_v2000_to_molecule(&v2000).unwrap();
        assert_eq!(parsed.atoms.len(), 3);
        assert_eq!(parsed.bonds.len(), 2);
        assert_eq!(parsed.compute_formula_properties().get_formula(), "C2H6O");

        // V3000 roundtrip + SDF multi-record roundtrip
        let v3000 = molecule_to_molfile_v3000(&mol, "EthanolV3000");
        let parsed_v3 = molfile_to_molecule(&v3000).unwrap();
        assert_eq!(parsed_v3.atoms.len(), 3);
        assert_eq!(parsed_v3.bonds.len(), 2);
        assert_eq!(parsed_v3.compute_formula_properties().get_formula(), "C2H6O");

        let mut props = BTreeMap::new();
        props.insert("ID".to_string(), "REC-01".to_string());
        let sdf = molecules_to_sdf(&[
            SdfRecord {
                title: "Ethanol".to_string(),
                molecule: mol.clone(),
                properties: props.clone(),
            },
            SdfRecord {
                title: "Ethanol2".to_string(),
                molecule: mol,
                properties: BTreeMap::new(),
            },
        ]);
        let recs = sdf_to_molecules(&sdf).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].properties.get("ID").map(String::as_str), Some("REC-01"));
    }

    #[test]
    fn test_smiles_cdx28_extended_and_cdxml_budget() {
        // 1. SMILES parse & roundtrip (benzene, acetic acid, ammonium chloride salt)
        for &(smi, expected_formula) in &[
            ("c1ccccc1", "C6H6"),
            ("CC(=O)O", "C2H4O2"),
            ("[NH4+].[Cl-]", "H4ClN"),
            ("CC#N", "C2H3N"),
            ("C12CC1C2", "C4H6"),
        ] {
            let m = smiles_to_molecule(smi).unwrap();
            assert_eq!(
                m.compute_formula_properties().get_formula(),
                expected_formula,
                "SMILES formula mismatch for {smi}"
            );
            let out_smi = molecule_to_smiles(&m);
            let reparsed = smiles_to_molecule(&out_smi).unwrap();
            assert_eq!(reparsed.atoms.len(), m.atoms.len());
            assert_eq!(reparsed.bonds.len(), m.bonds.len());
            assert_eq!(reparsed.compute_formula_properties().get_formula(), expected_formula);
        }

        // 2. 28-byte historical CDX header + >= 65536-byte extended property roundtrip
        let mut raw = Vec::new();
        raw.extend_from_slice(CDX_MAGIC_PREFIX);
        raw.extend_from_slice(&[0u8; 16]); // 28-byte header total
        raw.extend_from_slice(&CDX_OBJ_DOCUMENT.to_le_bytes());
        raw.extend_from_slice(&1u32.to_le_bytes());
        // Extended property payload (0xffff sentinel + u32 len = 70000)
        raw.extend_from_slice(&0x0700u16.to_le_bytes());
        raw.extend_from_slice(&0xffffu16.to_le_bytes());
        raw.extend_from_slice(&70_000u32.to_le_bytes());
        raw.extend(std::iter::repeat_n(b'A', 70_000));
        raw.extend_from_slice(&0u16.to_le_bytes()); // end Document
        raw.extend_from_slice(&0u16.to_le_bytes()); // file terminator

        let doc = CdxDocument::parse(&raw).unwrap();
        assert_eq!(doc.header_bytes.len(), 28);
        assert_eq!(doc.to_bytes(), raw);

        // 3. CDXML security & budget guards
        let xxe = r#"<!DOCTYPE foo [<!ENTITY xxe SYSTEM "file:///etc/passwd">]><CDXML/>"#;
        assert!(CdxmlDocument::parse(xxe).is_err());

        let tiny_budget = CdxmlParserBudget {
            max_bytes: 1024,
            max_elements: 2,
            max_depth: 16,
        };
        assert!(CdxmlDocument::parse_with_budget("<CDXML><page><fragment/></page></CDXML>", tiny_budget).is_err());
        let shallow_budget = CdxmlParserBudget {
            max_bytes: 1024,
            max_elements: 100,
            max_depth: 2,
        };
        assert!(CdxmlDocument::parse_with_budget("<CDXML><page><fragment/></page></CDXML>", shallow_budget).is_err());

        // 4. MDL RXN V2000 reaction roundtrip
        let r1 = smiles_to_molecule("C=C").unwrap();
        let p1 = smiles_to_molecule("CC").unwrap();
        let rxn = MdlRxnFile {
            title: "Hydrogenation".to_string(),
            reactants: vec![r1],
            products: vec![p1],
        };
        let rxn_txt = reaction_to_rxn_v2000(&rxn);
        let parsed_rxn = rxn_v2000_to_reaction(&rxn_txt).unwrap();
        assert_eq!(parsed_rxn.title, "Hydrogenation");
        assert_eq!(parsed_rxn.reactants.len(), 1);
        assert_eq!(parsed_rxn.products.len(), 1);
        assert_eq!(parsed_rxn.reactants[0].bonds[0].order, 2);
        assert_eq!(parsed_rxn.products[0].bonds[0].order, 1);

        // 5. Verify M RAD roundtrip & CDXML EnhancedStereo corpus extraction
        let mut rad_mol = ChemMolecule::new();
        let mut rad_c = ChemAtom::new(1, 6);
        rad_c.radical = 1; // Doublet (1e) in ChiralCraft binary enum
        rad_mol.add_atom(rad_c);
        let rad_v2k = molecule_to_molfile_v2000(&rad_mol, "MethylRadical");
        assert!(rad_v2k.contains("M  RAD  1   1   2"));
        let rad_back = molfile_to_molecule(&rad_v2k).unwrap();
        assert_eq!(rad_back.atoms[0].radical, 1);

        let enh_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../research/application/_7AE2F42284C3DE988F0C8A1221FFE197");
        if enh_path.exists() {
            let xml = fs::read_to_string(enh_path).unwrap();
            let cdxml = CdxmlDocument::parse(&xml).unwrap();
            let enh = cdxml.extract_enhanced_stereo();
            assert!(!enh.is_empty());
            let all_markers: Vec<_> = enh.values().flatten().copied().collect();
            assert!(all_markers.iter().any(|m| m.group_type == 4 && m.group_number == 1));
        }

        // 6. MDL RXN V3000 roundtrip + live ChiralCraft $RXN V3000 atom map & RXCTR check
        let v3_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../research/runtime-probes/p15/p15-mapped-rxn.rxn");
        if v3_path.exists() {
            let v3_txt = fs::read_to_string(v3_path).unwrap();
            let parsed_v3 = rxn_v3000_to_reaction(&v3_txt).unwrap();
            assert_eq!(parsed_v3.base.reactants.len(), 1);
            assert_eq!(parsed_v3.base.products.len(), 1);
            assert_eq!(parsed_v3.reactant_atom_maps, vec![vec![1, 2]]);
            assert_eq!(parsed_v3.product_atom_maps, vec![vec![1, 2, 0]]);
            assert_eq!(parsed_v3.reactant_bond_rxctrs, vec![vec![8]]);
            assert_eq!(parsed_v3.product_bond_rxctrs, vec![vec![8, 4]]);
            let reemitted = reaction_to_rxn_v3000(&parsed_v3);
            let reparsed = rxn_v3000_to_reaction(&reemitted).unwrap();
            assert_eq!(reparsed.reactant_atom_maps, parsed_v3.reactant_atom_maps);
            assert_eq!(reparsed.product_atom_maps, parsed_v3.product_atom_maps);
            assert_eq!(reparsed.reactant_bond_rxctrs, parsed_v3.reactant_bond_rxctrs);
            assert_eq!(reparsed.product_bond_rxctrs, parsed_v3.product_bond_rxctrs);
            assert_eq!(rxn_to_reaction(&reemitted).unwrap().products[0].atoms.len(), 3);
        }
    }
}
