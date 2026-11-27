//! Zero-Dependency Stdio MCP (Model Context Protocol) Server for ChiralCraft (`chiralcraft --mcp`)
//!
//! Exposes ChiralCraft's pure-Rust chemistry engine, NMR predictor, IUPAC Structure↔Name,
//! Editable SVG generator, PubChem/ChEMBL search, and file analyzer/converter over JSON-RPC 2.0 stdio
//! so AI coding and research agents can create, analyze, and render publication-ready chemistry directly.

use std::fs;
use std::io::{self, BufRead, Write};

use chiralcraft_chem::{
    compute_analysis_report, estimate_chemical_properties, molecule_to_name, predict_13c_nmr,
    predict_1h_nmr, NameDictionary,
};
use chiralcraft_doc::{
    clean_up_structure_preserving_stereo, DocumentSettings, InteractiveDocument,
};
use chiralcraft_formats::{molecule_to_molfile_v3000, molecule_to_smiles, smiles_to_molecule};

use crate::chem_db::{query_chembl_public, query_pubchem_public, DbSearchMode, JsonValue};
use crate::{load_document_any, save_document_any};

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn format_id(id: Option<&JsonValue>) -> String {
    match id {
        Some(JsonValue::Number(n)) => {
            if n.fract().abs() < 1e-9 {
                format!("{}", *n as i64)
            } else {
                format!("{n}")
            }
        }
        Some(JsonValue::String(s)) => format!("\"{}\"", json_escape(s)),
        _ => "null".to_string(),
    }
}

fn rpc_ok(id: Option<&JsonValue>, result_json: &str) -> String {
    format!(
        "{{\"jsonrpc\":\"2.0\",\"id\":{},\"result\":{result_json}}}",
        format_id(id)
    )
}

fn rpc_err(id: Option<&JsonValue>, code: i32, msg: &str) -> String {
    format!(
        "{{\"jsonrpc\":\"2.0\",\"id\":{},\"error\":{{\"code\":{code},\"message\":\"{}\"}}}}",
        format_id(id),
        json_escape(msg)
    )
}

fn tool_text_result(text: &str, is_error: bool) -> String {
    format!(
        "{{\"content\":[{{\"type\":\"text\",\"text\":\"{}\"}}],\"isError\":{}}}",
        json_escape(text),
        if is_error { "true" } else { "false" }
    )
}

fn resolve_stylesheet(name: Option<&str>) -> DocumentSettings {
    let key = name.unwrap_or("ACS").trim().to_ascii_lowercase();
    for (sname, st) in DocumentSettings::builtin_journal_stylesheets() {
        if sname.to_ascii_lowercase().contains(&key) {
            return st;
        }
    }
    DocumentSettings::acs_document_1996()
}

fn tools_list_json() -> &'static str {
    concat!(
        r#"{"tools":["#,
        r#"{"name":"analyze_smiles","description":"Analyze a SMILES string: formula, exact mass, MW, elemental analysis, LogP, tPSA, pKa, LogS, IUPAC name, and chemical warnings.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"predict_nmr","description":"Predict 1H-NMR (400 MHz) or 13C-NMR (100 MHz) chemical shifts and multiplicities from SMILES.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"nucleus":{"type":"string","enum":["1H","13C"]}},"required":["smiles"]}},"#,
        r#"{"name":"name_to_structure","description":"Convert an IUPAC or common chemical name into SMILES, MOL V3000, and an editable ChiralCraft SVG.","inputSchema":{"type":"object","properties":{"name":{"type":"string"},"stylesheet":{"type":"string"}},"required":["name"]}},"#,
        r#"{"name":"render_editable_svg","description":"Render a SMILES structure to a self-contained editable ChiralCraft SVG (compatible with WordCraft/DeckCraft) and optionally write it to output_path.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"stylesheet":{"type":"string"},"transparent":{"type":"boolean"},"dark_ink":{"type":"boolean"},"caption":{"type":"string"},"output_path":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"search_chem_db","description":"Search PubChem and/or ChEMBL by compound name, CAS, CID, ChEMBL ID, formula, or SMILES.","inputSchema":{"type":"object","properties":{"query":{"type":"string"},"source":{"type":"string","enum":["both","pubchem","chembl"]}},"required":["query"]}},"#,
        r#"{"name":"inspect_file","description":"Inspect a chemical or editable media file (.cdxml, .cdx, .mol, .sdf, .xyz, .cif, .pdb, POSCAR, .svg, .png, .pdf) and report molecules, formula, MW, and name.","inputSchema":{"type":"object","properties":{"path":{"type":"string"}},"required":["path"]}},"#,
        r#"{"name":"convert_file","description":"Convert between supported chemical/vector formats (.cdxml, .cdx, .mol, .sdf, .rxn, .xyz, .cif, .pdb, POSCAR, .svg, .png, .pdf, .eps).","inputSchema":{"type":"object","properties":{"input_path":{"type":"string"},"output_path":{"type":"string"}},"required":["input_path","output_path"]}},"#,
        r#"{"name":"diff_files","description":"Perform a chemistry-aware diff between two chemical or editable media files (.cdxml, .svg, .png, .pdf, .mol, .sdf).","inputSchema":{"type":"object","properties":{"left_path":{"type":"string"},"right_path":{"type":"string"}},"required":["left_path","right_path"]}},"#,
        r#"{"name":"medchem_scorecard","description":"Compute Lipinski, Veber, Ghose, CNS MPO, QED, Fsp3, SA score, and PAINS/Brenk alerts for a SMILES string.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"simulate_hrms","description":"Simulate HRMS isotopic envelope peaks and ESI+/ESI- adduct m/z table for a SMILES string.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"adduct":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"retrosynthesis_suggest","description":"Suggest strategic retrosynthetic disconnections, synthons, and forward reagents for a target SMILES.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"predict_spectroscopy","description":"Predict FT-IR bands (cm^-1), UV-Vis lambda_max & color, and 2D NMR (COSY/HSQC/HMBC) cross-peaks from SMILES.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"analyze_organometallic","description":"Analyze transition-metal complexes using the 18-electron rule & Green CBC (oxidation state, d^n, total valence e-, geometry).","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"export_latex","description":"Export a SMILES structure to compilable LaTeX \\chemfig{...} and TikZ vector code.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"generate_si_procedure","description":"Generate a publication-ready ACS/JOC Supporting Information experimental procedure + 1H/13C NMR, IR, and HRMS characterization block.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"name":{"type":"string"},"solvent":{"type":"string"},"yield_pct":{"type":"number"}},"required":["smiles"]}},"#,
        r#"{"name":"export_quantum_input","description":"Generate a ready-to-run ORCA (.inp), Gaussian (.gjf), or GFN2-xTB input deck (with explicit 3D hydrogens + UFF pre-opt) from SMILES.","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"program":{"type":"string","enum":["orca","gaussian","xtb"]},"method":{"type":"string"},"basis":{"type":"string"},"solvent":{"type":"string"}},"required":["smiles"]}},"#,
        r#"{"name":"parse_quantum_output","description":"Parse an ORCA .out, Gaussian .log, xTB g98.out/vibspectrum, or Molden file for SCF energy, HOMO/LUMO gap, dipole moment, imaginary modes, and harmonic IR frequencies + km/mol intensities.","inputSchema":{"type":"object","properties":{"path":{"type":"string"}},"required":["path"]}},"#,
        r#"{"name":"parse_jcamp_ir","description":"Parse an experimental FT-IR JCAMP-DX (.jdx/.dx, including ASDF/SQZ/DIF/DUP compression) or CSV/DPT file, pick prominent peaks, and optionally compute HQI/Pearson match against a SMILES.","inputSchema":{"type":"object","properties":{"path":{"type":"string"},"compare_smiles":{"type":"string"}},"required":["path"]}},"#,
        r#"{"name":"run_xtb_ir","description":"Run a local GFN2-xTB or GFN-FF harmonic IR calculation (--ohess) on a SMILES structure if xtb is installed on PATH, returning scaled frequencies and intensities (km/mol).","inputSchema":{"type":"object","properties":{"smiles":{"type":"string"},"method":{"type":"string","enum":["gfn2","gfnff"]},"solvent":{"type":"string"}},"required":["smiles"]}}"#,
        r#"]}"#
    )
}

/// Processes one JSON-RPC 2.0 MCP request line and returns `Some(response_line)` (or `None` for notifications).
pub fn handle_mcp_request(line: &str) -> Option<String> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return None;
    }
    let root = match JsonValue::parse(trimmed) {
        Ok(v) => v,
        Err(e) => return Some(rpc_err(None, -32700, &format!("Parse error: {e}"))),
    };
    let id = root.get("id");
    let method = root.get("method").and_then(|m| m.as_str()).unwrap_or("");

    match method {
        "notifications/initialized" | "notifications/cancelled" => None,
        "initialize" => Some(rpc_ok(
            id,
            concat!(
                r#"{"protocolVersion":"2025-06-18","capabilities":{"tools":{}},"serverInfo":{"name":"chiralcraft","version":""#,
                env!("CARGO_PKG_VERSION"),
                r#""}}"#
            ),
        )),
        "ping" => Some(rpc_ok(id, "{}")),
        "tools/list" => Some(rpc_ok(id, tools_list_json())),
        "tools/call" => {
            let params = root.get("params");
            let tname = params
                .and_then(|p| p.get("name"))
                .and_then(|n| n.as_str())
                .unwrap_or("");
            let args = params.and_then(|p| p.get("arguments"));
            let res = dispatch_tool(tname, args);
            Some(rpc_ok(id, &res))
        }
        other => Some(rpc_err(id, -32601, &format!("Method not found: {other}"))),
    }
}

fn dispatch_tool(name: &str, args: Option<&JsonValue>) -> String {
    let get_s = |k: &str| -> Option<&str> { args?.get(k)?.as_str() };
    let get_b = |k: &str, default: bool| -> bool {
        match args.and_then(|a| a.get(k)) {
            Some(JsonValue::Bool(b)) => *b,
            _ => default,
        }
    };
    let dict = NameDictionary::with_builtins();

    match name {
        "analyze_smiles" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let rep = compute_analysis_report(&mol);
                    let props = estimate_chemical_properties(&mol);
                    let iupac = molecule_to_name(&mol, smi, &dict);
                    let mut doc = InteractiveDocument::new_saved_seed("mcp", DocumentSettings::acs_document_1996());
                    doc.scene.molecules.push(mol);
                    let warns = doc.check_structure_with_preferences();
                    let summary = format!(
                        "Name: {iupac}\nFormula: {}\nExact Mass: {:.4}\nMol. Wt.: {:.2}\nElem. Anal.: {}\nLogP: {:.2} | tPSA: {:.2} | LogS: {:.3} | pKa: {:?}\nWarnings: {}",
                        rep.formula, rep.exact_mass, rep.mol_weight, rep.elem_analysis,
                        props.log_p, props.tpsa, props.log_s, props.pka, warns.len()
                    );
                    tool_text_result(&summary, false)
                }
                Err(e) => tool_text_result(&format!("SMILES error: {e}"), true),
            }
        }
        "predict_nmr" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let is_13c = get_s("nucleus").map(|n| n.eq_ignore_ascii_case("13c")).unwrap_or(false);
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let pred = if is_13c { predict_13c_nmr(&mol) } else { predict_1h_nmr(&mol) };
                    tool_text_result(&pred.report_text, false)
                }
                Err(e) => tool_text_result(&format!("SMILES error: {e}"), true),
            }
        }
        "name_to_structure" => {
            let Some(nm) = get_s("name") else {
                return tool_text_result("Missing required argument: name", true);
            };
            let st = resolve_stylesheet(get_s("stylesheet"));
            match dict.name_to_smiles(nm).and_then(|s| Ok((s.clone(), smiles_to_molecule(&s)?))) {
                Ok((smi, mut mol)) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], st.bond_length_pt);
                    let molfile = molecule_to_molfile_v3000(&mol, nm);
                    let mut doc = InteractiveDocument::new_saved_seed(nm, st);
                    doc.scene.molecules.push(mol);
                    let svg = doc.to_svg_transparent();
                    let text = format!("SMILES: {smi}\n\n--- MOL V3000 ---\n{molfile}\n--- EDITABLE SVG ---\n{svg}");
                    tool_text_result(&text, false)
                }
                Err(e) => tool_text_result(&format!("Name resolution failed: {e}"), true),
            }
        }
        "render_editable_svg" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let st = resolve_stylesheet(get_s("stylesheet"));
            let transparent = get_b("transparent", true);
            let dark_ink = get_b("dark_ink", false);
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], st.bond_length_pt);
                    let mut doc = InteractiveDocument::new_saved_seed("render.svg", st);
                    doc.scene.molecules.push(mol);
                    if let Some(cap_txt) = get_s("caption") {
                        doc.scene.captions.push(chiralcraft_doc::CaptionNode {
                            id: 900,
                            text: cap_txt.to_string(),
                            x: 0.0,
                            y: 42.0,
                            element: None,
                            generic_nickname: None,
                            warning: None,
                            justification: chiralcraft_doc::LabelJustification::Center,
                        });
                    }
                    let svg = doc.to_svg_styled(transparent, dark_ink, 1.0, 3.6);
                    if let Some(out_p) = get_s("output_path") {
                        if let Err(e) = fs::write(out_p, &svg) {
                            return tool_text_result(&format!("Write error for {out_p}: {e}"), true);
                        }
                        return tool_text_result(&format!("Wrote editable ChiralCraft SVG ({} bytes) to {out_p}", svg.len()), false);
                    }
                    tool_text_result(&svg, false)
                }
                Err(e) => tool_text_result(&format!("SMILES error: {e}"), true),
            }
        }
        "search_chem_db" => {
            let Some(q) = get_s("query") else {
                return tool_text_result("Missing required argument: query", true);
            };
            let src = get_s("source").unwrap_or("both").to_ascii_lowercase();
            let mut lines = Vec::new();
            if src == "both" || src == "pubchem" {
                if let Ok(hits) = query_pubchem_public(q, DbSearchMode::NameOrId, 85) {
                    for h in hits.into_iter().take(5) {
                        lines.push(format!("[PubChem] {} | {} | {} | MW {} | SMILES: {}", h.id, h.name, h.formula, h.molecular_weight, h.smiles));
                    }
                }
            }
            if src == "both" || src == "chembl" {
                if let Ok(hits) = query_chembl_public(q, DbSearchMode::NameOrId, 85) {
                    for h in hits.into_iter().take(5) {
                        let ph = h.max_phase.as_deref().unwrap_or("—");
                        lines.push(format!("[ChEMBL] {} | {} | {} | MW {} | Phase: {} | SMILES: {}", h.id, h.name, h.formula, h.molecular_weight, ph, h.smiles));
                    }
                }
            }
            if lines.is_empty() {
                tool_text_result("No matching compounds found.", false)
            } else {
                tool_text_result(&lines.join("\n"), false)
            }
        }
        "inspect_file" => {
            let Some(path) = get_s("path") else {
                return tool_text_result("Missing required argument: path", true);
            };
            match load_document_any(path) {
                Ok(doc) => {
                    let mut lines = vec![format!(
                        "File: {path} | Molecules: {} | Captions: {} | Arrows: {} | Stylesheet: {}",
                        doc.scene.molecules.len(),
                        doc.scene.captions.len(),
                        doc.scene.arrows.len(),
                        doc.settings.stylesheet_name
                    )];
                    for (i, m) in doc.scene.molecules.iter().enumerate() {
                        let smi = molecule_to_smiles(m);
                        let rep = compute_analysis_report(m);
                        let nm = molecule_to_name(m, &smi, &dict);
                        lines.push(format!("  #{}: {} | {} | MW {:.2} | SMILES: {}", i + 1, nm, rep.formula, rep.mol_weight, smi));
                    }
                    tool_text_result(&lines.join("\n"), false)
                }
                Err(e) => tool_text_result(&format!("Inspect failed: {e}"), true),
            }
        }
        "convert_file" => {
            let (Some(inp), Some(outp)) = (get_s("input_path"), get_s("output_path")) else {
                return tool_text_result("Missing input_path or output_path", true);
            };
            match load_document_any(inp).and_then(|d| save_document_any(&d, outp)) {
                Ok(()) => tool_text_result(&format!("Converted {inp} -> {outp}"), false),
                Err(e) => tool_text_result(&format!("Convert failed: {e}"), true),
            }
        }
        "diff_files" => {
            let (Some(lp), Some(rp)) = (get_s("left_path"), get_s("right_path")) else {
                return tool_text_result("Missing left_path or right_path", true);
            };
            match (load_document_any(lp), load_document_any(rp)) {
                (Ok(ld), Ok(rd)) => {
                    let diff = chiralcraft_doc::compare_documents(&ld, &rd);
                    tool_text_result(&diff.format_text_report(), false)
                }
                (Err(e), _) => tool_text_result(&format!("Failed to load {lp}: {e}"), true),
                (_, Err(e)) => tool_text_result(&format!("Failed to load {rp}: {e}"), true),
            }
        }
        "medchem_scorecard" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mol) => {
                    let mc = chiralcraft_chem::compute_medchem_report(&mol);
                    let alerts_txt = if mc.alerts.is_empty() {
                        "None".to_string()
                    } else {
                        mc.alerts.iter().map(|a| format!("[{}] {}", a.category, a.name)).collect::<Vec<_>>().join("; ")
                    };
                    let txt = format!(
                        "MW: {:.2} | LogP: {:.2} | tPSA: {:.1} | HBD: {} | HBA: {} | RotB: {} | Fsp3: {:.2}\nLipinski: {} ({} viol) | Veber: {} | Ghose: {}\nCNS MPO: {:.2}/6.0 | QED: {:.2} | SA Score: {:.2}/10\nPAINS/Brenk Alerts ({}): {alerts_txt}",
                        mc.mw, mc.log_p, mc.tpsa, mc.hbd, mc.hba, mc.rotatable_bonds, mc.fsp3,
                        if mc.lipinski_pass { "PASS" } else { "FAIL" }, mc.lipinski_violations,
                        if mc.veber_pass { "PASS" } else { "FAIL" },
                        if mc.ghose_pass { "PASS" } else { "FAIL" },
                        mc.cns_mpo, mc.qed, mc.sa_score, mc.alerts.len()
                    );
                    tool_text_result(&txt, false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "simulate_hrms" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let adduct = get_s("adduct").unwrap_or("[M+H]+");
            match smiles_to_molecule(smi) {
                Ok(mol) => {
                    let sim = chiralcraft_chem::simulate_hrms(&mol, adduct);
                    let mut lines = vec![format!("Formula: {} | Neutral Exact Mass: {:.4} | Adduct: {}", sim.formula, sim.neutral_exact_mass, sim.selected_adduct)];
                    lines.push("Isotope Peaks (m/z : rel%):".to_string());
                    for pk in sim.peaks.iter().take(8) {
                        lines.push(format!("  {:.4} : {:.1}%", pk.mz, pk.rel_intensity));
                    }
                    lines.push("Adducts:".to_string());
                    for ad in sim.adducts.iter().take(8) {
                        lines.push(format!("  {} = {:.4}", ad.name, ad.mz));
                    }
                    tool_text_result(&lines.join("\n"), false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "retrosynthesis_suggest" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let hits = chiralcraft_doc::suggest_disconnections(&mol, 14.40);
                    if hits.is_empty() {
                        tool_text_result("No acyclic strategic disconnections suggested.", false)
                    } else {
                        let lines: Vec<String> = hits
                            .iter()
                            .enumerate()
                            .map(|(i, h)| format!("{}. {} (score {}) -> {} | Reagents: {}", i + 1, h.reaction_name, h.score, h.synthons_summary, h.reagents))
                            .collect();
                        tool_text_result(&lines.join("\n"), false)
                    }
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "predict_spectroscopy" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mol) => {
                    let ir = chiralcraft_chem::predict_ir_spectrum(&mol);
                    let uv = chiralcraft_chem::predict_uv_vis(&mol);
                    let hsqc = chiralcraft_chem::predict_2d_nmr(&mol, chiralcraft_chem::Nmr2dKind::Hsqc);
                    let hmbc = chiralcraft_chem::predict_2d_nmr(&mol, chiralcraft_chem::Nmr2dKind::Hmbc);
                    let txt = format!(
                        "{}\nUV-Vis: λmax ≈ {:.0} nm (ε ≈ {:.0}), {} — {}\n2D NMR: {} HSQC cross-peaks, {} HMBC cross-peaks",
                        ir.report_text, uv.lambda_max_nm, uv.epsilon_est, uv.transition, uv.perceived_color_name,
                        hsqc.peaks.len(), hmbc.peaks.len()
                    );
                    tool_text_result(&txt, false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "analyze_organometallic" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mol) => {
                    let reps = chiralcraft_chem::analyze_organometallic_complexes(&mol);
                    if reps.is_empty() {
                        tool_text_result("No d-block transition metal centers found.", false)
                    } else {
                        let lines: Vec<String> = reps
                            .iter()
                            .map(|r| format!("{} (Group {}): {}e-, ox {:+}, d^{}, ML{}X{}, {} — {}", r.metal_symbol, r.group_number, r.total_valence_electrons, r.oxidation_state, r.d_electrons, r.l_ligands, r.x_ligands, r.predicted_geometry, r.rule_status))
                            .collect();
                        tool_text_result(&lines.join("\n"), false)
                    }
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "export_latex" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let cf = chiralcraft_doc::molecule_to_chemfig(&mol);
                    let mut doc = InteractiveDocument::new_saved_seed("tex", DocumentSettings::acs_document_1996());
                    doc.scene.molecules.push(mol);
                    let tikz = chiralcraft_doc::document_to_tikz(&doc);
                    tool_text_result(&format!("{cf}\n\n{tikz}"), false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "generate_si_procedure" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let cname = get_s("name").unwrap_or("Target Compound");
            let solv = get_s("solvent").unwrap_or("THF");
            let yld = args
                .and_then(|a| a.get("yield_pct"))
                .and_then(|v| match v { JsonValue::Number(n) => Some(*n), _ => None })
                .unwrap_or(90.0);
            match smiles_to_molecule(smi) {
                Ok(mol) => {
                    let si = chiralcraft_doc::generate_si_experimental_block(&mol, cname, &[], solv, 5.0, "80 °C", "12 h", yld);
                    tool_text_result(&si, false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "export_quantum_input" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let mut spec = chiralcraft_formats::QuantumJobSpec::default();
            spec.program = match get_s("program").unwrap_or("orca").to_ascii_lowercase().as_str() {
                "gaussian" | "gjf" => chiralcraft_formats::QuantumProgram::Gaussian,
                "xtb" => chiralcraft_formats::QuantumProgram::Xtb,
                _ => chiralcraft_formats::QuantumProgram::Orca,
            };
            if let Some(m) = get_s("method") { spec.method = m.to_string(); }
            if let Some(b) = get_s("basis") { spec.basis = b.to_string(); }
            if let Some(s) = get_s("solvent") { spec.solvent = s.to_string(); }
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let lifted = chiralcraft_doc::lift_molecule_2d_to_3d(&mol);
                    let deck = chiralcraft_formats::export_quantum_input(&lifted, &spec);
                    tool_text_result(&deck, false)
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        "parse_quantum_output" => {
            let Some(path) = get_s("path") else {
                return tool_text_result("Missing required argument: path", true);
            };
            match fs::read_to_string(path)
                .map_err(|e| e.to_string())
                .and_then(|s| chiralcraft_formats::parse_quantum_output(&s))
            {
                Ok(rep) => {
                    let mut top_modes: Vec<&chiralcraft_formats::QuantumIrMode> = rep
                        .ir_modes
                        .iter()
                        .filter(|m| !m.is_imaginary)
                        .collect();
                    top_modes.sort_by(|a, b| {
                        b.intensity_km_mol
                            .partial_cmp(&a.intensity_km_mol)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let mode_lines: Vec<String> = top_modes
                        .iter()
                        .take(8)
                        .map(|m| format!("  Mode #{}: {:.1} cm^-1 (Int {:.1} km/mol)", m.mode_index, m.wavenumber_cm1, m.intensity_km_mol))
                        .collect();
                    let txt = format!(
                        "Program: {} | Converged: {}\nSCF Energy (Eh): {:?}\nHOMO (eV): {:?} | LUMO (eV): {:?} | Gap (eV): {:?}\nDipole (Debye): {:?}\nIR Modes: {} (Imaginary: {})\n{}",
                        rep.program,
                        rep.converged,
                        rep.energy_hartree,
                        rep.homo_ev,
                        rep.lumo_ev,
                        rep.homo_lumo_gap_ev,
                        rep.dipole_debye,
                        rep.ir_modes.len().max(rep.ir_frequencies_cm1.len()),
                        rep.imaginary_mode_count,
                        mode_lines.join("\n")
                    );
                    tool_text_result(&txt, false)
                }
                Err(e) => tool_text_result(&format!("Quantum parse error: {e}"), true),
            }
        }
        "parse_jcamp_ir" => {
            let Some(path) = get_s("path") else {
                return tool_text_result("Missing required argument: path", true);
            };
            match fs::read_to_string(path)
                .map_err(|e| e.to_string())
                .and_then(|s| chiralcraft_formats::parse_experimental_ir_file(path, &s))
            {
                Ok(spec) => {
                    let abs = spec.ir_absorbance_curve();
                    let peaks = chiralcraft_chem::pick_ir_peaks(&abs, 0.12, 10);
                    let pk_str: Vec<String> = peaks.iter().map(|p| format!("{:.0} cm^-1", p.wavenumber_cm1)).collect();
                    let mut txt = format!(
                        "Title: {} | IR: {} | Points: {}\nTop Experimental Peaks: {}",
                        spec.title,
                        spec.is_infrared(),
                        spec.raw_points.len(),
                        pk_str.join(", ")
                    );
                    if let Some(smi) = get_s("compare_smiles") {
                        if let Ok(mol) = smiles_to_molecule(smi) {
                            let ir = chiralcraft_chem::predict_ir_spectrum(&mol);
                            let params = chiralcraft_chem::IrRenderParams::default();
                            let calc_abs = ir.simulate_absorbance(&params, 400);
                            let cmp = chiralcraft_chem::compare_ir_spectra(&calc_abs, &abs, 400.0, 4000.0);
                            txt.push_str(&format!(
                                "\nComparison vs '{}' (Rule IR): HQI = {:.1}% | Pearson r = {:.3}",
                                smi, cmp.hqi_pct, cmp.pearson_r
                            ));
                        }
                    }
                    tool_text_result(&txt, false)
                }
                Err(e) => tool_text_result(&format!("Experimental IR parse error: {e}"), true),
            }
        }
        "run_xtb_ir" => {
            let Some(smi) = get_s("smiles") else {
                return tool_text_result("Missing required argument: smiles", true);
            };
            let bins = crate::quantum_runner::detect_qc_binaries("", "", "");
            let Some(xtb_bin) = bins.xtb_path else {
                return tool_text_result("xtb binary not found on PATH. Install via `conda install -c conda-forge xtb` or export an input deck with `export_quantum_input`.", true);
            };
            let backend = match get_s("method").unwrap_or("gfn2").to_ascii_lowercase().as_str() {
                "gfnff" | "gfn-ff" => crate::quantum_runner::QcRunBackend::XtbGfnFf,
                _ => crate::quantum_runner::QcRunBackend::XtbGfn2,
            };
            let mut spec = chiralcraft_formats::QuantumJobSpec::default();
            if let Some(s) = get_s("solvent") {
                spec.solvent = s.to_string();
            }
            match smiles_to_molecule(smi) {
                Ok(mut mol) => {
                    clean_up_structure_preserving_stereo(&mut mol, &[], 14.40);
                    let lifted = chiralcraft_doc::lift_molecule_2d_to_3d(&mol);
                    let cancel = std::sync::atomic::AtomicBool::new(false);
                    match crate::quantum_runner::run_qc_job_blocking(
                        backend,
                        &lifted,
                        format!("{smi}|{backend:?}"),
                        &spec,
                        &xtb_bin,
                        &cancel,
                    ) {
                        Ok(out) => {
                            let scale = out.ir_prediction.default_scale_factor;
                            let mut lines = vec![format!(
                                "{} completed in {:.1}s | Modes: {} | Imaginary: {}\n{}",
                                out.ir_prediction.source.label(),
                                out.elapsed_secs,
                                out.ir_prediction.bands.len(),
                                out.ir_prediction.imaginary_count,
                                out.ir_prediction.format_report(scale)
                            )];
                            for b in out.ir_prediction.bands.iter().take(12) {
                                lines.push(format!(
                                    "  {:.1} cm^-1 (Int {:.1} km/mol, {}) — {}",
                                    b.scaled_wavenumber(scale),
                                    b.intensity_km_mol.unwrap_or(0.0),
                                    b.strength,
                                    b.assignment
                                ));
                            }
                            tool_text_result(&lines.join("\n"), false)
                        }
                        Err(e) => tool_text_result(&format!("xTB run failed: {e}"), true),
                    }
                }
                Err(e) => tool_text_result(&format!("Invalid SMILES: {e}"), true),
            }
        }
        other => tool_text_result(&format!("Unknown tool: {other}"), true),
    }
}

/// Runs the line-delimited JSON-RPC 2.0 MCP server on stdin/stdout (`chiralcraft --mcp`).
pub fn run_mcp_stdio_server() {
    let stdin = io::stdin();
    let mut stdout = io::stdout();
    for line_res in stdin.lock().lines() {
        let Ok(line) = line_res else { break };
        if let Some(reply) = handle_mcp_request(&line) {
            let _ = writeln!(stdout, "{reply}");
            let _ = stdout.flush();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mcp_protocol_initialize_list_and_tools() {
        let init = handle_mcp_request(r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}"#).unwrap();
        assert!(!init.contains('\n'), "MCP NDJSON responses must be single-line");
        assert!(init.contains("\"chiralcraft\"") && init.contains("protocolVersion"));
        assert!(JsonValue::parse(&init).is_ok());

        let list = handle_mcp_request(r#"{"jsonrpc":"2.0","id":2,"method":"tools/list"}"#).unwrap();
        assert!(!list.contains('\n'), "tools/list must be single-line NDJSON");
        let parsed = JsonValue::parse(&list).expect("tools/list must be valid JSON");
        let tools = parsed.get("result").and_then(|r| r.get("tools")).and_then(|t| t.as_array()).unwrap();
        assert_eq!(tools.len(), 19, "All 19 MCP tools must be exposed");

        let analyze = handle_mcp_request(
            r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"analyze_smiles","arguments":{"smiles":"CCO"}}}"#,
        ).unwrap();
        assert!(!analyze.contains('\n') && JsonValue::parse(&analyze).is_ok());
        assert!(analyze.contains("C2H6O") && analyze.contains("ethanol"));

        let svg = handle_mcp_request(
            r#"{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"render_editable_svg","arguments":{"smiles":"c1ccccc1","stylesheet":"Nature","dark_ink":true}}}"#,
        ).unwrap();
        assert!(!svg.contains('\n') && JsonValue::parse(&svg).is_ok());
        assert!(svg.contains("chiralcraft:document") && svg.contains("#F2F5F9"));

        // Verify parse_jcamp_ir end-to-end with HQI comparison
        let tmp_jdx = std::env::temp_dir().join(format!("chiralcraft-mcp-ir-{}.jdx", std::process::id()));
        let sample_jdx = "##TITLE=Test Ketone IR\n##DATA TYPE=INFRARED SPECTRUM\n##XUNITS=1/CM\n##YUNITS=TRANSMITTANCE\n##FIRSTX=4000\n##LASTX=400\n##NPOINTS=5\n##XYDATA=(X++(Y..Y))\n4000 0.98 0.95 0.15 0.92 0.97\n##END=\n";
        fs::write(&tmp_jdx, sample_jdx).unwrap();
        let jdx_req = format!(
            r#"{{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{{"name":"parse_jcamp_ir","arguments":{{"path":"{}","compare_smiles":"CC(=O)C"}}}}}}"#,
            json_escape(&tmp_jdx.to_string_lossy())
        );
        let jdx_resp = handle_mcp_request(&jdx_req).unwrap();
        let _ = fs::remove_file(&tmp_jdx);
        assert!(!jdx_resp.contains('\n') && JsonValue::parse(&jdx_resp).is_ok());
        assert!(jdx_resp.contains("Test Ketone IR") && jdx_resp.contains("HQI ="));
    }
}
