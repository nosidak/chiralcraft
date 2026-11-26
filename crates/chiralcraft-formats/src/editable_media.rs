//! Lossless Editable PNG (`tEXt` metadata chunk) & Native Vector PDF 1.4 Export/Import.
//!
//! Enables true round-trip editing not just for SVG, but also for raster PNG figures
//! (embedding the lossless CDXML chemical graph inside a standard PNG `tEXt` chunk
//! `ChiralCraftDocument`) and publication-grade vector PDFs (pure-Rust PDF 1.4 with
//! embedded lossless `%CHIRALCRAFT_CDXML_B64:` payload).

use chiralcraft_chem::{get_element, BondDisplay, CcVal, ChemMolecule};

const PNG_SIG: [u8; 8] = [137, 80, 78, 71, 13, 10, 26, 10];
const PNG_TEXT_KEYWORD: &[u8] = b"ChiralCraftDocument";
const PDF_B64_PREFIX: &str = "%CHIRALCRAFT_CDXML_B64:";

fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFF_u32;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            if (crc & 1) != 0 {
                crc = (crc >> 1) ^ 0xEDB8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}

/// Standard RFC 4648 Base64 encoder (zero dependencies).
pub fn base64_encode(input: &[u8]) -> String {
    const TBL: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(((input.len() + 2) / 3) * 4);
    for chunk in input.chunks(3) {
        let b0 = chunk[0] as u32;
        let b1 = if chunk.len() > 1 { chunk[1] as u32 } else { 0 };
        let b2 = if chunk.len() > 2 { chunk[2] as u32 } else { 0 };
        let n = (b0 << 16) | (b1 << 8) | b2;
        out.push(TBL[((n >> 18) & 63) as usize] as char);
        out.push(TBL[((n >> 12) & 63) as usize] as char);
        out.push(if chunk.len() > 1 { TBL[((n >> 6) & 63) as usize] as char } else { '=' });
        out.push(if chunk.len() > 2 { TBL[(n & 63) as usize] as char } else { '=' });
    }
    out
}

/// Standard RFC 4648 Base64 decoder (ignores ASCII whitespace).
pub fn base64_decode(input: &str) -> Result<Vec<u8>, String> {
    let val = |c: u8| -> Option<u32> {
        match c {
            b'A'..=b'Z' => Some((c - b'A') as u32),
            b'a'..=b'z' => Some((c - b'a' + 26) as u32),
            b'0'..=b'9' => Some((c - b'0' + 52) as u32),
            b'+' => Some(62),
            b'/' => Some(63),
            _ => None,
        }
    };
    let clean: Vec<u8> = input
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .collect();
    if clean.len() % 4 != 0 {
        return Err("Invalid base64 length".to_string());
    }
    let mut out = Vec::with_capacity((clean.len() / 4) * 3);
    for chunk in clean.chunks(4) {
        let v0 = val(chunk[0]).ok_or_else(|| "Bad base64 char".to_string())?;
        let v1 = val(chunk[1]).ok_or_else(|| "Bad base64 char".to_string())?;
        let pad2 = chunk[2] == b'=';
        let pad3 = chunk[3] == b'=';
        let v2 = if pad2 { 0 } else { val(chunk[2]).ok_or_else(|| "Bad base64 char".to_string())? };
        let v3 = if pad3 { 0 } else { val(chunk[3]).ok_or_else(|| "Bad base64 char".to_string())? };
        let n = (v0 << 18) | (v1 << 12) | (v2 << 6) | v3;
        out.push(((n >> 16) & 0xFF) as u8);
        if !pad2 {
            out.push(((n >> 8) & 0xFF) as u8);
        }
        if !pad3 {
            out.push((n & 0xFF) as u8);
        }
    }
    Ok(out)
}

/// Injects a standard PNG `tEXt` chunk (`ChiralCraftDocument\0<base64(cdxml)>`) right before `IEND`.
/// The output remains a 100% compliant PNG readable by every browser, image viewer, and office suite.
pub fn embed_chiralcraft_cdxml_in_png(png_bytes: &[u8], cdxml_xml: &str) -> Vec<u8> {
    if png_bytes.len() < 20 || &png_bytes[..8] != PNG_SIG {
        return png_bytes.to_vec();
    }
    let encoded = base64_encode(cdxml_xml.as_bytes());
    let mut payload = Vec::with_capacity(PNG_TEXT_KEYWORD.len() + 1 + encoded.len());
    payload.extend_from_slice(PNG_TEXT_KEYWORD);
    payload.push(0);
    payload.extend_from_slice(encoded.as_bytes());

    let mut text_chunk = Vec::with_capacity(12 + payload.len());
    text_chunk.extend_from_slice(&(payload.len() as u32).to_be_bytes());
    text_chunk.extend_from_slice(b"tEXt");
    text_chunk.extend_from_slice(&payload);
    let c = crc32(&text_chunk[4..]);
    text_chunk.extend_from_slice(&c.to_be_bytes());

    // Walk chunks and splice before IEND (replacing any previous ChiralCraftDocument chunk)
    let mut out = Vec::with_capacity(png_bytes.len() + text_chunk.len());
    out.extend_from_slice(&PNG_SIG);
    let mut pos = 8usize;
    let mut inserted = false;
    while pos + 12 <= png_bytes.len() {
        let len = u32::from_be_bytes([
            png_bytes[pos],
            png_bytes[pos + 1],
            png_bytes[pos + 2],
            png_bytes[pos + 3],
        ]) as usize;
        let ctype = &png_bytes[pos + 4..pos + 8];
        let total = 12 + len;
        if pos + total > png_bytes.len() {
            break;
        }
        if ctype == b"tEXt" && png_bytes[pos + 8..pos + 8 + len].starts_with(PNG_TEXT_KEYWORD) {
            pos += total;
            continue;
        }
        if ctype == b"IEND" && !inserted {
            out.extend_from_slice(&text_chunk);
            inserted = true;
        }
        out.extend_from_slice(&png_bytes[pos..pos + total]);
        pos += total;
    }
    if !inserted {
        out.extend_from_slice(&text_chunk);
    }
    out
}

/// Extracts an embedded ChiralCraft CDXML document from a PNG `tEXt` chunk, if present.
pub fn extract_chiralcraft_cdxml_from_png(png_bytes: &[u8]) -> Option<String> {
    if png_bytes.len() < 20 || &png_bytes[..8] != PNG_SIG {
        return None;
    }
    let mut pos = 8usize;
    while pos + 12 <= png_bytes.len() {
        let len = u32::from_be_bytes([
            png_bytes[pos],
            png_bytes[pos + 1],
            png_bytes[pos + 2],
            png_bytes[pos + 3],
        ]) as usize;
        let ctype = &png_bytes[pos + 4..pos + 8];
        let total = 12 + len;
        if pos + total > png_bytes.len() {
            break;
        }
        if ctype == b"tEXt" {
            let data = &png_bytes[pos + 8..pos + 8 + len];
            if let Some(nul) = data.iter().position(|&b| b == 0) {
                if &data[..nul] == PNG_TEXT_KEYWORD {
                    let b64 = std::str::from_utf8(&data[nul + 1..]).ok()?;
                    let decoded = base64_decode(b64).ok()?;
                    return String::from_utf8(decoded).ok();
                }
            }
        }
        pos += total;
    }
    None
}

/// Lightweight caption descriptor for vector PDF rendering.
#[derive(Debug, Clone)]
pub struct PdfCaptionItem {
    pub x: f64,
    pub y: f64,
    pub text: String,
}

/// Lightweight reaction arrow descriptor for vector PDF rendering.
#[derive(Debug, Clone)]
pub struct PdfArrowItem {
    pub tail: (f64, f64),
    pub head: (f64, f64),
}

fn escape_pdf_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            c if c.is_ascii() && !c.is_ascii_control() => out.push(c),
            'β' => out.push_str("beta"),
            'α' => out.push_str("alpha"),
            'δ' => out.push_str("d"),
            '–' | '—' => out.push('-'),
            '₀' => out.push('0'),
            '₁' => out.push('1'),
            '₂' => out.push('2'),
            '₃' => out.push('3'),
            '₄' => out.push('4'),
            '₅' => out.push('5'),
            '₆' => out.push('6'),
            '₇' => out.push('7'),
            '₈' => out.push('8'),
            '₉' => out.push('9'),
            '⁺' => out.push('+'),
            '⁻' => out.push('-'),
            _ => out.push(' '),
        }
    }
    out
}

/// Renders molecules, arrows, and captions to a pure-Rust publication-grade vector **PDF 1.4**
/// containing an embedded lossless ChiralCraft CDXML payload for round-trip reopening!
pub fn scene_to_editable_pdf(
    molecules: &[ChemMolecule],
    arrows: &[PdfArrowItem],
    captions: &[PdfCaptionItem],
    cdxml_xml: &str,
    title: &str,
    dark_ink: bool,
    padding_pt: f64,
) -> Vec<u8> {
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
    for m in molecules {
        for a in &m.atoms {
            if a.visible {
                expand(a.x, a.y);
            }
        }
    }
    for ar in arrows {
        expand(ar.tail.0, ar.tail.1);
        expand(ar.head.0, ar.head.1);
    }
    for c in captions {
        expand(c.x, c.y);
    }
    if !min_x.is_finite() {
        min_x = 0.0;
        min_y = 0.0;
        max_x = 240.0;
        max_y = 140.0;
    }
    let pad = padding_pt.max(12.0);
    let page_w = ((max_x - min_x) + 2.0 * pad).max(100.0);
    let page_h = ((max_y - min_y) + 2.0 * pad).max(80.0);
    let tx = |x: f64| (x - min_x) + pad;
    let ty = |y: f64| page_h - ((y - min_y) + pad); // PDF origin is bottom-left!

    let mut content = String::new();
    if dark_ink {
        content.push_str("0.118 0.133 0.165 rg\n");
        content.push_str(&format!("0 0 {page_w:.2} {page_h:.2} re f\n"));
        content.push_str("0.949 0.961 0.976 RG 0.949 0.961 0.976 rg\n");
    } else {
        content.push_str("0 0 0 RG 0 0 0 rg\n");
    }
    content.push_str("0.85 w 1 J 1 j\n");

    for mol in molecules {
        let val = CcVal::compute(mol);
        for b in &mol.bonds {
            if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
                continue;
            }
            let a1 = &mol.atoms[b.begin_atom];
            let a2 = &mol.atoms[b.end_atom];
            let (x1, y1) = (tx(a1.x), ty(a1.y));
            let (x2, y2) = (tx(a2.x), ty(a2.y));
            let dx = x2 - x1;
            let dy = y2 - y1;
            let len = dx.hypot(dy).max(1e-3);
            let nx = -dy / len;
            let ny = dx / len;
            match b.display {
                BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                    let (sx, sy, ex, ey) = if b.display == BondDisplay::WedgeBegin {
                        (x1, y1, x2, y2)
                    } else {
                        (x2, y2, x1, y1)
                    };
                    let hw = 2.0;
                    content.push_str(&format!(
                        "{sx:.2} {sy:.2} m {:.2} {:.2} l {:.2} {:.2} l h f\n",
                        ex + nx * hw, ey + ny * hw, ex - nx * hw, ey - ny * hw
                    ));
                }
                BondDisplay::Bold => {
                    content.push_str(&format!("2.1 w {x1:.2} {y1:.2} m {x2:.2} {y2:.2} l S 0.85 w\n"));
                }
                _ => {
                    if b.order == 2 {
                        let off = 1.4;
                        content.push_str(&format!(
                            "{:.2} {:.2} m {:.2} {:.2} l S\n{:.2} {:.2} m {:.2} {:.2} l S\n",
                            x1 + nx * off, y1 + ny * off, x2 + nx * off, y2 + ny * off,
                            x1 - nx * off, y1 - ny * off, x2 - nx * off, y2 - ny * off,
                        ));
                    } else if b.order == 3 {
                        let off = 2.1;
                        content.push_str(&format!(
                            "{x1:.2} {y1:.2} m {x2:.2} {y2:.2} l S\n{:.2} {:.2} m {:.2} {:.2} l S\n{:.2} {:.2} m {:.2} {:.2} l S\n",
                            x1 + nx * off, y1 + ny * off, x2 + nx * off, y2 + ny * off,
                            x1 - nx * off, y1 - ny * off, x2 - nx * off, y2 - ny * off,
                        ));
                    } else {
                        content.push_str(&format!("{x1:.2} {y1:.2} m {x2:.2} {y2:.2} l S\n"));
                    }
                }
            }
        }
        for (ai, a) in mol.atoms.iter().enumerate() {
            if !a.visible {
                continue;
            }
            let lbl = if let Some(l) = &a.label {
                Some(l.clone())
            } else if a.z != 6 || a.charge != 0 {
                let sym = get_element(a.z).map(|e| e.symbol).unwrap_or("?");
                let h = val.implicit_h_counts.get(ai).copied().unwrap_or(0);
                let mut s = sym.to_string();
                if h == 1 { s.push('H'); } else if h > 1 { s.push_str(&format!("H{h}")); }
                if a.charge > 0 { s.push('+'); } else if a.charge < 0 { s.push('-'); }
                Some(s)
            } else {
                None
            };
            if let Some(txt) = lbl {
                let px = tx(a.x);
                let py = ty(a.y);
                let bw = (txt.len() as f64 * 5.2).max(8.0);
                // Knock out underlying bond junction with background fill
                if dark_ink {
                    content.push_str(&format!("0.118 0.133 0.165 rg {:.2} {:.2} {:.2} 9.5 re f 0.949 0.961 0.976 rg\n", px - bw * 0.5, py - 4.5, bw));
                } else {
                    content.push_str(&format!("1 1 1 rg {:.2} {:.2} {:.2} 9.5 re f 0 0 0 rg\n", px - bw * 0.5, py - 4.5, bw));
                }
                content.push_str(&format!(
                    "BT /F1 9.5 Tf {:.2} {:.2} Td ({}) Tj ET\n",
                    px - bw * 0.42,
                    py - 3.2,
                    escape_pdf_str(&txt)
                ));
            }
        }
    }

    for ar in arrows {
        let (x1, y1) = (tx(ar.tail.0), ty(ar.tail.1));
        let (x2, y2) = (tx(ar.head.0), ty(ar.head.1));
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = dx.hypot(dy).max(1e-3);
        let ux = dx / len;
        let uy = dy / len;
        let nx = -uy;
        let ny = ux;
        content.push_str(&format!("{x1:.2} {y1:.2} m {x2:.2} {y2:.2} l S\n"));
        content.push_str(&format!(
            "{x2:.2} {y2:.2} m {:.2} {:.2} l {:.2} {:.2} l h f\n",
            x2 - ux * 5.5 + nx * 2.2, y2 - uy * 5.5 + ny * 2.2,
            x2 - ux * 5.5 - nx * 2.2, y2 - uy * 5.5 - ny * 2.2,
        ));
    }
    for cap in captions {
        let px = tx(cap.x);
        let py = ty(cap.y);
        content.push_str(&format!(
            "BT /F1 9.5 Tf {:.2} {:.2} Td ({}) Tj ET\n",
            px,
            py - 3.0,
            escape_pdf_str(&cap.text)
        ));
    }

    let cdxml_b64 = base64_encode(cdxml_xml.as_bytes());
    let mut pdf = Vec::new();
    pdf.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
    let mut offsets = vec![0usize];

    let push_obj = |pdf: &mut Vec<u8>, offsets: &mut Vec<usize>, id: usize, body: &str| {
        offsets.push(pdf.len());
        pdf.extend_from_slice(format!("{id} 0 obj\n{body}\nendobj\n").as_bytes());
    };

    push_obj(&mut pdf, &mut offsets, 1, "<< /Type /Catalog /Pages 2 0 R >>");
    push_obj(&mut pdf, &mut offsets, 2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    push_obj(
        &mut pdf,
        &mut offsets,
        3,
        &format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {page_w:.2} {page_h:.2}] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
        ),
    );
    push_obj(
        &mut pdf,
        &mut offsets,
        4,
        "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    );
    push_obj(
        &mut pdf,
        &mut offsets,
        5,
        &format!("<< /Length {} >>\nstream\n{}endstream", content.len(), content),
    );
    push_obj(
        &mut pdf,
        &mut offsets,
        6,
        &format!(
            "<< /Title ({}) /Creator (ChiralCraft 1.0) /Producer (ChiralCraft Vector PDF Engine) >>",
            escape_pdf_str(title)
        ),
    );

    let xref_pos = pdf.len();
    pdf.extend_from_slice(format!("xref\n0 {}\n0000000000 65535 f \n", offsets.len()).as_bytes());
    for &off in offsets.iter().skip(1) {
        pdf.extend_from_slice(format!("{off:010} 00000 n \n").as_bytes());
    }
    pdf.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R /Info 6 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n{PDF_B64_PREFIX}{cdxml_b64}\n",
            offsets.len()
        )
        .as_bytes(),
    );
    pdf
}

/// Extracts an embedded lossless ChiralCraft CDXML document from an exported PDF, if present.
pub fn extract_chiralcraft_cdxml_from_pdf(pdf_bytes: &[u8]) -> Option<String> {
    let text = String::from_utf8_lossy(pdf_bytes);
    let idx = text.rfind(PDF_B64_PREFIX)?;
    let rest = &text[idx + PDF_B64_PREFIX.len()..];
    let b64 = rest.lines().next()?.trim();
    let decoded = base64_decode(b64).ok()?;
    String::from_utf8(decoded).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{molecule_to_png, smiles_to_molecule};

    #[test]
    fn test_editable_png_and_pdf_lossless_roundtrip() {
        let mol = smiles_to_molecule("CCO").expect("ethanol");
        let raw_png = molecule_to_png(&mol, 2.0);
        let sample_cdxml = "<CDXML><page><fragment id=\"1\"/></page></CDXML>";
        let ed_png = embed_chiralcraft_cdxml_in_png(&raw_png, sample_cdxml);
        assert_eq!(&ed_png[..8], &PNG_SIG);
        assert!(crate::decode_png_rgba(&ed_png).is_some(), "Editable PNG must still decode cleanly");
        assert_eq!(
            extract_chiralcraft_cdxml_from_png(&ed_png).as_deref(),
            Some(sample_cdxml)
        );

        let ed_pdf = scene_to_editable_pdf(
            &[mol],
            &[PdfArrowItem { tail: (10.0, 10.0), head: (50.0, 10.0) }],
            &[PdfCaptionItem { x: 20.0, y: 30.0, text: "Ethanol".into() }],
            sample_cdxml,
            "Ethanol",
            false,
            18.0,
        );
        assert!(ed_pdf.starts_with(b"%PDF-1.4"));
        assert!(String::from_utf8_lossy(&ed_pdf).contains("%%EOF"));
        assert_eq!(
            extract_chiralcraft_cdxml_from_pdf(&ed_pdf).as_deref(),
            Some(sample_cdxml)
        );
    }
}
