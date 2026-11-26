//! Modular File I/O, Atomic Round-Trip Writer, and Linux Wayland/X11 Editable SVG Clipboard Bridge.
//!
//! Handles 22+ chemical, crystallographic, quantum, and lossless editable media formats:
//! `.cdxml`, `.cdx`, `.cds`, `.mol`, `.sdf`, `.rxn`, `.ct`, `.cml`, `.tgf`, `.skc`,
//! `.xyz`, `.cif`, `.pdb`, `POSCAR`/`CONTCAR`, `.inp`/`.gjf`/`.out`/`.log`,
//! editable `.svg`, editable `.png`, editable `.pdf`, `.eps`, `.emf`, `.sln`, `.inchi`.

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

use chiralcraft_doc::{
    clean_up_structure_preserving_stereo, svg_contains_chiralcraft_document, CopyBackground,
    DocumentSettings, ExportOptions, InteractiveDocument,
};
use chiralcraft_formats::{
    cif_to_molecule_with_cell, cml_to_molecules, ct_to_molecule, inchi_to_inchikey,
    molecule_to_cml, molecule_to_ct, molecule_to_emf, molecule_to_eps, molecule_to_inchi,
    molecule_to_molfile_v3000, molecule_to_skc, molecule_to_sln, molecule_to_tgf, molecules_to_sdf,
    molfile_to_molecule, molfile_v3000_to_molecule, pdb_to_molecule, poscar_to_molecule_with_cell,
    quantum_file_to_molecule, reaction_to_rxn_v2000, rxn_to_reaction, sdf_to_molecules,
    skc_to_molecules, smiles_to_molecule, tgf_to_molecules, xyz_to_molecule, CdxmlDocument,
    CdxDocument, CrystalUnitCell, MdlRxnFile, SdfRecord,
};

pub fn atomic_write_text(path: &str, content: &str) -> Result<(), String> {
    let tmp = format!("{path}.tmp.{}", std::process::id());
    fs::write(&tmp, content).map_err(|e| e.to_string())?;
    fs::rename(&tmp, path)
        .or_else(|_| {
            let _ = fs::remove_file(&tmp);
            fs::write(path, content)
        })
        .map_err(|e| e.to_string())
}

/// Pushes an Editable SVG (`image/svg+xml`) onto the Linux Wayland (`wl-copy`) or X11 (`xclip`) clipboard
/// so `Ctrl+C` in ChiralCraft pastes directly as a vector graphic into WordCraft / DeckCraft / LibreOffice.
pub fn copy_svg_to_os_clipboard(svg: &str) -> Result<&'static str, String> {
    if cfg!(test) {
        return Err("Test environment".to_string());
    }
    let has_wayland = env::var_os("WAYLAND_DISPLAY").is_some();
    let has_x11 = env::var_os("DISPLAY").is_some();
    if has_wayland {
        if let Ok(mut child) = Command::new("wl-copy")
            .args(["--type", "image/svg+xml"])
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
        {
            if let Some(mut stdin) = child.stdin.take() {
                let _ = stdin.write_all(svg.as_bytes());
            }
            let _ = child.wait();
            return Ok("Wayland image/svg+xml");
        }
    }
    if has_x11 {
        if let Ok(mut child) = Command::new("xclip")
            .args(["-selection", "clipboard", "-t", "image/svg+xml", "-i"])
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
        {
            if let Some(mut stdin) = child.stdin.take() {
                let _ = stdin.write_all(svg.as_bytes());
            }
            let _ = child.wait();
            return Ok("X11 image/svg+xml");
        }
    }
    Err("No Wayland/X11 SVG clipboard tool available".to_string())
}

/// Reads an Editable ChiralCraft SVG from the OS clipboard (`image/svg+xml` or SVG text).
pub fn paste_svg_from_os_clipboard() -> Option<String> {
    if cfg!(test) {
        return None;
    }
    if env::var_os("WAYLAND_DISPLAY").is_some() {
        for args in [vec!["--type", "image/svg+xml"], vec![]] {
            if let Ok(out) = Command::new("wl-paste").args(&args).output() {
                if out.status.success() {
                    let s = String::from_utf8_lossy(&out.stdout).to_string();
                    if svg_contains_chiralcraft_document(&s) {
                        return Some(s);
                    }
                }
            }
        }
    }
    if env::var_os("DISPLAY").is_some() {
        for args in [
            vec!["-selection", "clipboard", "-t", "image/svg+xml", "-o"],
            vec!["-selection", "clipboard", "-o"],
        ] {
            if let Ok(out) = Command::new("xclip").args(&args).output() {
                if out.status.success() {
                    let s = String::from_utf8_lossy(&out.stdout).to_string();
                    if svg_contains_chiralcraft_document(&s) {
                        return Some(s);
                    }
                }
            }
        }
    }
    None
}

pub fn load_document_any_with_cell(
    path: &str,
) -> Result<(InteractiveDocument, Option<CrystalUnitCell>), String> {
    let low = path.to_ascii_lowercase();
    let fname = std::path::Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if low.ends_with(".svg") {
        let s = fs::read_to_string(path).map_err(|e| e.to_string())?;
        return Ok((InteractiveDocument::from_svg(&s)?, None));
    } else if low.ends_with(".png") {
        let bytes = fs::read(path).map_err(|e| e.to_string())?;
        return Ok((InteractiveDocument::from_png(&bytes)?, None));
    } else if low.ends_with(".pdf") {
        let bytes = fs::read(path).map_err(|e| e.to_string())?;
        return Ok((InteractiveDocument::from_pdf(&bytes)?, None));
    } else if low.ends_with(".cdx") || low.ends_with(".cds") {
        let bytes = fs::read(path).map_err(|e| e.to_string())?;
        let cdxml = CdxDocument::parse(&bytes)
            .map_err(|e| format!("{e:?}"))?
            .to_cdxml();
        Ok((InteractiveDocument::from_cdxml(&cdxml), None))
    } else if low.ends_with(".cdxml") {
        let s = fs::read_to_string(path).map_err(|e| e.to_string())?;
        let cdxml = CdxmlDocument::parse(&s).map_err(|e| format!("{e:?}"))?;
        Ok((InteractiveDocument::from_cdxml(&cdxml), None))
    } else {
        let mut doc =
            InteractiveDocument::new_saved_seed(path, DocumentSettings::acs_document_1996());
        let mut cell_opt = None;
        let mols = if low.ends_with(".sdf") || low.ends_with(".sd") {
            sdf_to_molecules(&fs::read_to_string(path).map_err(|e| e.to_string())?)?
                .into_iter()
                .map(|r| r.molecule)
                .collect()
        } else if low.ends_with(".rxn") {
            let rxn = rxn_to_reaction(&fs::read_to_string(path).map_err(|e| e.to_string())?)?;
            let mut all = rxn.reactants;
            all.extend(rxn.products);
            all
        } else if low.ends_with(".mol") {
            let s = fs::read_to_string(path).map_err(|e| e.to_string())?;
            vec![if s.contains("V3000") {
                molfile_v3000_to_molecule(&s)?
            } else {
                molfile_to_molecule(&s)?
            }]
        } else if low.ends_with(".ct") {
            vec![ct_to_molecule(&fs::read_to_string(path).map_err(|e| e.to_string())?)?]
        } else if low.ends_with(".cml") {
            cml_to_molecules(&fs::read_to_string(path).map_err(|e| e.to_string())?)?
        } else if low.ends_with(".tgf") {
            tgf_to_molecules(&fs::read_to_string(path).map_err(|e| e.to_string())?)?
        } else if low.ends_with(".skc") {
            skc_to_molecules(&fs::read(path).map_err(|e| e.to_string())?)?
        } else if low.ends_with(".xyz") {
            vec![xyz_to_molecule(&fs::read_to_string(path).map_err(|e| e.to_string())?)?]
        } else if low.ends_with(".cif") {
            let (m, cell) =
                cif_to_molecule_with_cell(&fs::read_to_string(path).map_err(|e| e.to_string())?)?;
            cell_opt = Some(cell);
            vec![m]
        } else if low.ends_with(".pdb") {
            vec![pdb_to_molecule(&fs::read_to_string(path).map_err(|e| e.to_string())?)?]
        } else if fname == "poscar"
            || fname == "contcar"
            || low.ends_with(".poscar")
            || low.ends_with(".vasp")
        {
            let (m, cell) = poscar_to_molecule_with_cell(
                &fs::read_to_string(path).map_err(|e| e.to_string())?,
            )?;
            cell_opt = Some(cell);
            vec![m]
        } else if low.ends_with(".out")
            || low.ends_with(".log")
            || low.ends_with(".inp")
            || low.ends_with(".gjf")
        {
            vec![quantum_file_to_molecule(
                &fs::read_to_string(path).map_err(|e| e.to_string())?,
            )?]
        } else if low.ends_with(".smi") || low.ends_with(".smiles") {
            let txt = fs::read_to_string(path).map_err(|e| e.to_string())?;
            let mut out_mols = Vec::new();
            for line in txt.lines() {
                let t = line.trim();
                if t.is_empty() || t.starts_with('#') {
                    continue;
                }
                if let Some(tok) = t.split_whitespace().next() {
                    if let Ok(mut m) = smiles_to_molecule(tok) {
                        clean_up_structure_preserving_stereo(&mut m, &[], 14.40);
                        out_mols.push(m);
                    }
                }
            }
            if out_mols.is_empty() {
                return Err("No valid SMILES entries found in file".to_string());
            }
            out_mols
        } else {
            return Err("Unsupported input format".to_string());
        };
        for m in mols {
            doc.scene.molecules.push(m);
        }
        Ok((doc, cell_opt))
    }
}

pub fn load_document_any(path: &str) -> Result<InteractiveDocument, String> {
    load_document_any_with_cell(path).map(|(d, _)| d)
}

/// Loads either a file path on disk (any supported ChiralCraft format) or, if `input` is not
/// a file path, parses it directly as a SMILES string and lays out clean 2D coordinates.
pub fn load_document_or_smiles(input: &str) -> Result<InteractiveDocument, String> {
    let p = std::path::Path::new(input);
    if p.exists() {
        return load_document_any(input);
    }
    let looks_like_file = input.contains('/')
        || [
            ".cdxml", ".cdx", ".cds", ".mol", ".sdf", ".sd", ".rxn", ".ct", ".cml", ".tgf",
            ".skc", ".xyz", ".cif", ".pdb", ".svg", ".png", ".pdf", ".out", ".log", ".inp",
            ".gjf", ".smi", ".smiles", ".poscar", ".vasp",
        ]
        .iter()
        .any(|ext| input.to_ascii_lowercase().ends_with(ext));
    if looks_like_file {
        return load_document_any(input);
    }
    let mut mol = smiles_to_molecule(input)?;
    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
    let mut doc = InteractiveDocument::new_saved_seed(input, DocumentSettings::acs_document_1996());
    doc.scene.molecules.push(mol);
    Ok(doc)
}

pub fn save_document_any_with_opts(
    doc: &InteractiveDocument,
    out_path: &str,
    opts: ExportOptions,
    dark_ink: bool,
) -> Result<(), String> {
    let low = out_path.to_ascii_lowercase();
    if low.ends_with(".svg") {
        let transparent = matches!(opts.copy_bg, CopyBackground::Transparent);
        return atomic_write_text(
            out_path,
            &doc.to_svg_styled(transparent, dark_ink, 0.45, 3.6),
        );
    }
    if low.ends_with(".pdf") {
        return fs::write(out_path, doc.to_editable_pdf(dark_ink, opts.svg_padding_pt))
            .map_err(|e| e.to_string());
    }
    if low.ends_with(".png") {
        return fs::write(out_path, doc.to_editable_png(opts.png_dpi)?).map_err(|e| e.to_string());
    }
    if low.ends_with(".cdx") || low.ends_with(".cds") {
        return fs::write(out_path, doc.to_cdxml().to_cdx().to_bytes()).map_err(|e| e.to_string());
    }
    if low.ends_with(".cdxml") {
        return fs::write(out_path, doc.to_cdxml().to_xml_string()).map_err(|e| e.to_string());
    }
    if low.ends_with(".sdf") || low.ends_with(".sd") {
        let recs: Vec<SdfRecord> = doc
            .scene
            .molecules
            .iter()
            .enumerate()
            .map(|(i, m)| SdfRecord {
                title: format!("{}_{}", doc.name, i + 1),
                molecule: m.clone(),
                properties: BTreeMap::new(),
            })
            .collect();
        return fs::write(out_path, molecules_to_sdf(&recs)).map_err(|e| e.to_string());
    }
    if low.ends_with(".rxn") {
        let split_x = doc
            .scene
            .arrows
            .first()
            .map(|a| 0.5 * (a.tail.0 + a.head.0))
            .unwrap_or(240.0);
        let mut reactants = Vec::new();
        let mut products = Vec::new();
        for m in &doc.scene.molecules {
            let cx = if m.atoms.is_empty() {
                0.0
            } else {
                m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64)
            };
            if cx <= split_x {
                reactants.push(m.clone());
            } else {
                products.push(m.clone());
            }
        }
        let rxn = MdlRxnFile {
            title: doc.name.clone(),
            reactants,
            products,
        };
        return fs::write(out_path, reaction_to_rxn_v2000(&rxn)).map_err(|e| e.to_string());
    }
    let mol = doc
        .scene
        .molecules
        .first()
        .ok_or_else(|| "No molecule on canvas".to_string())?;
    if low.ends_with(".mol") {
        fs::write(out_path, molecule_to_molfile_v3000(mol, &doc.name)).map_err(|e| e.to_string())
    } else if low.ends_with(".ct") {
        fs::write(out_path, molecule_to_ct(mol, &doc.name)).map_err(|e| e.to_string())
    } else if low.ends_with(".tgf") {
        fs::write(out_path, molecule_to_tgf(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".cml") {
        fs::write(out_path, molecule_to_cml(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".eps") {
        fs::write(out_path, molecule_to_eps(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".emf") {
        fs::write(out_path, molecule_to_emf(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".skc") {
        fs::write(out_path, molecule_to_skc(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".sln") {
        fs::write(out_path, molecule_to_sln(mol)).map_err(|e| e.to_string())
    } else if low.ends_with(".inchi") {
        let inchi = molecule_to_inchi(mol);
        let key = inchi_to_inchikey(&inchi);
        fs::write(out_path, format!("{inchi}\n{key}\n")).map_err(|e| e.to_string())
    } else {
        fs::write(out_path, doc.to_cdxml().to_xml_string()).map_err(|e| e.to_string())
    }
}

pub fn save_document_any(doc: &InteractiveDocument, out_path: &str) -> Result<(), String> {
    save_document_any_with_opts(doc, out_path, ExportOptions::default(), false)
}
