//! Pure-Rust IUPAC Structure ↔ Name Engine for ChiralCraft parity
//! (`Structure > Convert Structure to Name`, `Convert Name to Structure`,
//!  `Add Structure to Dictionary...`, `Edit > Paste Special > Name as Structure`).
//!
//! Supports:
//! - 70+ canonical trivial, solvent, amino acid, heterocycle & pharmaceutical names
//! - Extensible runtime user dictionary (`NameDictionary`)
//! - Systematic IUPAC Name → SMILES / Structure parser (locants, multipliers, substituents,
//!   alkanes/alkenes/alkynes, cycloalkanes, substituted benzenes/pyridines, alcohols,
//!   aldehydes, ketones, carboxylic acids, esters, amides, amines, nitriles, thiols)
//! - Systematic Structure → IUPAC Name generator with functional group priority,
//!   ring/chain parent selection, lowest-locant numbering, and `(R)`/`(S)`/`(E)`/`(Z)` stereo.

use std::collections::{BTreeMap, HashMap};

use crate::{
    perceive_stereochemistry, BondDisplay, CcAromaticityInfo, CcVal, ChemAtom, ChemBond,
    ChemMolecule,
};

#[derive(Debug, Clone, Default)]
pub struct NameDictionary {
    /// Lowercase name -> SMILES
    pub name_to_smiles: BTreeMap<String, String>,
    /// Canonical formula + connectivity fingerprint -> preferred display name
    pub user_entries: Vec<(String, String)>,
}

impl NameDictionary {
    pub fn with_builtins() -> Self {
        let mut d = Self::default();
        for &(name, smi) in BUILTIN_NAMES {
            d.name_to_smiles.insert(name.to_lowercase(), smi.to_string());
        }
        d
    }

    pub fn add_entry(&mut self, name: &str, smiles: &str) {
        let clean = name.trim().to_string();
        if !clean.is_empty() && !smiles.trim().is_empty() {
            self.name_to_smiles.insert(clean.to_lowercase(), smiles.trim().to_string());
            self.user_entries.push((clean, smiles.trim().to_string()));
        }
    }

    /// Resolves any chemical name (trivial, pharmaceutical, user-defined, OR systematic IUPAC)
    /// into a valid SMILES string ready for 2D structure generation & cleanup.
    pub fn name_to_smiles(&self, input: &str) -> Result<String, String> {
        let raw = input.trim();
        if raw.is_empty() {
            return Err("Empty chemical name".to_string());
        }
        let low = raw.to_lowercase();
        if let Some(smi) = self.name_to_smiles.get(&low) {
            return Ok(smi.clone());
        }
        // Try systematic IUPAC parser!
        parse_systematic_iupac_to_smiles(&low)
            .ok_or_else(|| format!("Unrecognized chemical name: \"{raw}\""))
    }
}

/// Curated ChiralCraft built-in trivial, reagent, heterocycle & drug name table.
pub const BUILTIN_NAMES: &[(&str, &str)] = &[
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("acetylsalicylic acid", "CC(=O)Oc1ccccc1C(=O)O"),
    ("salicylic acid", "O=C(O)c1ccccc1O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("anisole", "COc1ccccc1"),
    ("benzoic acid", "O=C(O)c1ccccc1"),
    ("benzaldehyde", "O=Cc1ccccc1"),
    ("acetophenone", "CC(=O)c1ccccc1"),
    ("benzophenone", "O=C(c1ccccc1)c1ccccc1"),
    ("styrene", "C=Cc1ccccc1"),
    ("benzonitrile", "N#Cc1ccccc1"),
    ("nitrobenzene", "O=[N+]([O-])c1ccccc1"),
    ("chlorobenzene", "Clc1ccccc1"),
    ("bromobenzene", "Brc1ccccc1"),
    ("fluorobenzene", "Fc1ccccc1"),
    ("iodobenzene", "Ic1ccccc1"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
    ("norbornane", "C1CC2CCC1C2"),
    ("cubane", "C12C3C4C1C5C2C3C45"),
    ("decalin", "C1CCC2CCCCC2C1"),
    ("acetanilide", "CC(=O)Nc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("pyrimidine", "c1cncnc1"),
    ("pyrazine", "c1cnccn1"),
    ("pyrrole", "c1cc[nH]c1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("imidazole", "c1c[nH]cn1"),
    ("oxazole", "c1cocn1"),
    ("thiazole", "c1cscn1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("isoquinoline", "c1ccc2cnccc2c1"),
    ("purine", "c1ncc2[nH]cnc2n1"),
    ("piperidine", "C1CCNCC1"),
    ("piperazine", "C1CNCCN1"),
    ("morpholine", "C1COCCN1"),
    ("pyrrolidine", "C1CCNC1"),
    ("tetrahydrofuran", "C1CCOC1"),
    ("thf", "C1CCOC1"),
    ("dioxane", "C1COCCO1"),
    ("water", "O"),
    ("ammonia", "N"),
    ("methane", "C"),
    ("ethane", "CC"),
    ("propane", "CCC"),
    ("butane", "CCCC"),
    ("pentane", "CCCCC"),
    ("hexane", "CCCCCC"),
    ("heptane", "CCCCCCC"),
    ("octane", "CCCCCCCC"),
    ("ethylene", "C=C"),
    ("acetylene", "C#C"),
    ("methanol", "CO"),
    ("ethanol", "CCO"),
    ("isopropanol", "CC(O)C"),
    ("glycerol", "OCC(O)CO"),
    ("formaldehyde", "C=O"),
    ("acetaldehyde", "CC=O"),
    ("acetone", "CC(=O)C"),
    ("formic acid", "O=CO"),
    ("acetic acid", "CC(=O)O"),
    ("oxalic acid", "O=C(O)C(=O)O"),
    ("citric acid", "O=C(O)CC(O)(C(=O)O)CC(=O)O"),
    ("ethyl acetate", "CCOC(=O)C"),
    ("diethyl ether", "CCOCC"),
    ("acetonitrile", "CC#N"),
    ("chloroform", "ClC(Cl)Cl"),
    ("dichloromethane", "ClCCl"),
    ("dcm", "ClCCl"),
    ("dmf", "CN(C)C=O"),
    ("dimethylformamide", "CN(C)C=O"),
    ("dmso", "CS(=O)C"),
    ("dimethyl sulfoxide", "CS(=O)C"),
    ("urea", "NC(=O)N"),
    ("triethylamine", "CCN(CC)CC"),
    ("glycine", "NCC(=O)O"),
    ("alanine", "CC(N)C(=O)O"),
    ("phenylalanine", "NC(Cc1ccccc1)C(=O)O"),
    ("dopamine", "NCCc1ccc(O)c(O)c1"),
    ("serotonin", "NCCc1c[nH]c2ccc(O)cc12"),
    ("nicotine", "CN1CCCC1c1cccnc1"),
    ("penicillin g", "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"),
    ("glucose", "OCC1OC(O)C(O)C(O)C1O"),
];

const CHAIN_ROOTS: &[(&str, usize)] = &[
    ("meth", 1),
    ("eth", 2),
    ("prop", 3),
    ("but", 4),
    ("pent", 5),
    ("hex", 6),
    ("hept", 7),
    ("oct", 8),
    ("non", 9),
    ("dec", 10),
    ("undec", 11),
    ("dodec", 12),
];

fn subst_to_smiles(name: &str) -> Option<&'static str> {
    match name {
        "fluoro" => Some("F"),
        "chloro" => Some("Cl"),
        "bromo" => Some("Br"),
        "iodo" => Some("I"),
        "hydroxy" | "hydroxyl" => Some("O"),
        "amino" => Some("N"),
        "mercapto" | "sulfanyl" => Some("S"),
        "nitro" => Some("[N+](=O)[O-]"),
        "cyano" => Some("C#N"),
        "methyl" => Some("C"),
        "ethyl" => Some("CC"),
        "propyl" => Some("CCC"),
        "isopropyl" => Some("C(C)C"),
        "butyl" => Some("CCCC"),
        "tert-butyl" | "t-butyl" | "tertbutyl" => Some("C(C)(C)C"),
        "methoxy" => Some("OC"),
        "ethoxy" => Some("OCC"),
        "phenyl" => Some("c1ccccc1"),
        "acetyl" => Some("C(=O)C"),
        "formyl" => Some("C=O"),
        "carboxy" => Some("C(=O)O"),
        "trifluoromethyl" => Some("C(F)(F)F"),
        "oxo" => Some("=O"),
        _ => None,
    }
}

fn strip_multiplier(s: &str) -> &str {
    for p in ["tetra", "tri", "di"] {
        if let Some(rest) = s.strip_prefix(p) {
            if subst_to_smiles(rest).is_some() {
                return rest;
            }
        }
    }
    s
}

/// Parses systematic IUPAC names like `2-chlorobutane`, `4-nitrobenzoic acid`, `2,4-dimethylphenol`,
/// `3-pentanone`, `hexanoic acid`, `2-methylpropan-1-ol`, `cyclohexanol`, `2-bromopyridine`,
/// `but-2-ene`, `prop-1-yne`, `ethyl propanoate` into SMILES!
fn parse_systematic_iupac_to_smiles(raw: &str) -> Option<String> {
    let mut s = raw.trim().to_string();
    // Strip leading stereo descriptors (R)-, (S)-, (E)-, (Z)-, cis-, trans-
    for pref in ["(r)-", "(s)-", "(e)-", "(z)-", "(2r)-", "(2s)-", "cis-", "trans-"] {
        if let Some(rest) = s.strip_prefix(pref) {
            s = rest.to_string();
        }
    }
    // Expand ortho/meta/para prefixes -> 2-/3-/4-
    if let Some(rest) = s.strip_prefix("ortho-").or_else(|| s.strip_prefix("o-")) {
        s = format!("2-{rest}");
    } else if let Some(rest) = s.strip_prefix("meta-").or_else(|| s.strip_prefix("m-")) {
        s = format!("3-{rest}");
    } else if let Some(rest) = s.strip_prefix("para-").or_else(|| s.strip_prefix("p-")) {
        s = format!("4-{rest}");
    }

    // 1. Esters: "<alkyl> <alkanoate>" (e.g. "methyl benzoate", "ethyl propanoate")
    if let Some((alk, acid_part)) = s.split_once(' ') {
        if let Some(alk_root) = alk.strip_suffix("yl") {
            if let Some(&(_, alk_n)) = CHAIN_ROOTS.iter().find(|&&(r, _)| r == alk_root) {
                let alk_smi = "C".repeat(alk_n);
                if acid_part == "benzoate" {
                    return Some(format!("{alk_smi}OC(=O)c1ccccc1"));
                }
                if acid_part == "formate" || acid_part == "methanoate" {
                    return Some(format!("{alk_smi}OC=O"));
                }
                if acid_part == "acetate" || acid_part == "ethanoate" {
                    return Some(format!("{alk_smi}OC(=O)C"));
                }
                if let Some(ac_root) = acid_part.strip_suffix("anoate") {
                    if let Some(&(_, ac_n)) = CHAIN_ROOTS.iter().find(|&&(r, _)| r == ac_root) {
                        let tail = "C".repeat(ac_n.saturating_sub(1));
                        return Some(format!("{alk_smi}OC(=O){tail}"));
                    }
                }
            }
        }
    }

    // 2. Substituted aromatic parents (`benzene`, `phenol`, `aniline`, `toluene`, `benzoic acid`, `benzaldehyde`, `pyridine`)
    let arom_parents: &[(&str, Option<&str>, bool)] = &[
        ("benzoic acid", Some("C(=O)O"), false),
        ("benzaldehyde", Some("C=O"), false),
        ("benzonitrile", Some("C#N"), false),
        ("phenol", Some("O"), false),
        ("aniline", Some("N"), false),
        ("toluene", Some("C"), false),
        ("anisole", Some("OC"), false),
        ("styrene", Some("C=C"), false),
        ("benzene", None, false),
        ("pyridine", None, true),
    ];
    for &(pname, pos1_sub, is_pyr) in arom_parents {
        if let Some(prefix) = s.strip_suffix(pname) {
            let Some(subs) = parse_substituent_prefix_checked(prefix) else {
                continue;
            };
            let mut ring_subs: [Vec<String>; 6] = Default::default();
            if let Some(p1) = pos1_sub {
                ring_subs[0].push(p1.to_string());
            }
            let mut next_ring = 2usize;
            for (loc, mut smi) in subs {
                if smi.contains('1') && next_ring <= 9 {
                    smi = smi.replace('1', &next_ring.to_string());
                    next_ring += 1;
                }
                let idx = loc.saturating_sub(1).min(5);
                ring_subs[idx].push(smi);
            }
            let mut out = String::new();
            for i in 0..6 {
                if i == 0 {
                    out.push_str(if is_pyr { "n1" } else { "c1" });
                } else {
                    out.push('c');
                }
                for sub in &ring_subs[i] {
                    out.push('(');
                    out.push_str(sub);
                    out.push(')');
                }
            }
            out.push('1');
            return Some(out);
        }
    }

    // 3. Cycloalkanes & substituted cycloalkane derivatives (`cyclohexanol`, `2-methylcyclopentanone`, `cyclobutane`)
    if let Some(c_idx) = s.find("cyclo") {
        let prefix = &s[..c_idx];
        let rest = &s[c_idx + 5..];
        for &(root, n) in CHAIN_ROOTS {
            if let Some(after_root) = rest.strip_prefix(root) {
                if n < 3 {
                    continue;
                }
                let Some(subs) = parse_substituent_prefix_checked(prefix) else {
                    continue;
                };
                let mut ring_subs: Vec<Vec<String>> = vec![Vec::new(); n];
                let mut dbl_pos: Option<usize> = None;
                match after_root {
                    "ane" => {}
                    "ene" => dbl_pos = Some(0),
                    "anol" => ring_subs[0].push("O".to_string()),
                    "anone" => ring_subs[0].push("=O".to_string()),
                    "anamine" => ring_subs[0].push("N".to_string()),
                    "anethiol" => ring_subs[0].push("S".to_string()),
                    "anecarboxylic acid" => ring_subs[0].push("C(=O)O".to_string()),
                    _ => continue,
                }
                let mut next_ring = 2usize;
                for (loc, mut smi) in subs {
                    if smi.contains('1') && next_ring <= 9 {
                        smi = smi.replace('1', &next_ring.to_string());
                        next_ring += 1;
                    }
                    let idx = loc.saturating_sub(1).min(n - 1);
                    ring_subs[idx].push(smi);
                }
                let mut out = String::new();
                for i in 0..n {
                    if i == 0 {
                        out.push_str("C1");
                    } else {
                        if dbl_pos == Some(i - 1) {
                            out.push('=');
                        }
                        out.push('C');
                    }
                    for sub in &ring_subs[i] {
                        if sub.starts_with('=') {
                            out.push('(');
                            out.push_str(sub);
                            out.push(')');
                        } else {
                            out.push('(');
                            out.push_str(sub);
                            out.push(')');
                        }
                    }
                }
                out.push('1');
                return Some(out);
            }
        }
    }

    // 4. Acyclic systematic IUPAC chains: [substituents][root][an/en/yn][suffix]
    for &(root, n) in CHAIN_ROOTS.iter().rev() {
        if let Some(rpos) = s.rfind(root) {
            let prefix = &s[..rpos];
            let suffix = &s[rpos + root.len()..];
            let Some((mut bond_orders, mut fn_group)) = parse_acyclic_suffix(suffix, n) else {
                continue;
            };
            // Support leading parent locant before root (e.g. `3-pentanone`, `2-butanol`, `2-butene`)
            let mut clean_prefix = prefix.trim().trim_end_matches('-');
            let mut parent_loc_override: Option<usize> = None;
            if let Some(last_tok) = clean_prefix.rsplit('-').next() {
                if let Ok(loc) = last_tok.parse::<usize>() {
                    parent_loc_override = Some(loc);
                    clean_prefix = clean_prefix.strip_suffix(last_tok).unwrap_or("").trim_end_matches('-');
                }
            }
            if let Some(loc) = parent_loc_override {
                if let Some((_, ref grp)) = fn_group {
                    fn_group = Some((loc, grp.clone()));
                } else if let Some(first_b) = bond_orders.first_mut() {
                    first_b.0 = loc.saturating_sub(1);
                }
            }
            let Some(subs) = parse_substituent_prefix_checked(clean_prefix) else {
                continue;
            };
            let mut atom_subs: Vec<Vec<String>> = vec![Vec::new(); n];
            let mut orders = vec![1i32; n.saturating_sub(1)];
            for (bpos, ord) in bond_orders {
                if bpos < orders.len() {
                    orders[bpos] = ord;
                }
            }
            if let Some((loc, grp)) = fn_group {
                let idx = loc.saturating_sub(1).min(n - 1);
                atom_subs[idx].push(grp);
            }
            for (loc, smi) in subs {
                let idx = loc.saturating_sub(1).min(n - 1);
                atom_subs[idx].push(smi);
            }
            let mut out = String::new();
            for i in 0..n {
                if i > 0 {
                    match orders[i - 1] {
                        2 => out.push('='),
                        3 => out.push('#'),
                        _ => {}
                    }
                }
                out.push('C');
                for sub in &atom_subs[i] {
                    out.push('(');
                    out.push_str(sub);
                    out.push(')');
                }
            }
            return Some(out);
        }
    }

    None
}

fn parse_acyclic_suffix(suffix: &str, chain_len: usize) -> Option<(Vec<(usize, i32)>, Option<(usize, String)>)> {
    let s = suffix.trim();
    match s {
        "ane" => Some((vec![], None)),
        "ene" => Some((vec![(0, 2)], None)),
        "yne" => Some((vec![(0, 3)], None)),
        "anol" => Some((vec![], Some((1, "O".to_string())))),
        "anal" => Some((vec![], Some((1, "=O".to_string())))),
        "anone" => Some((vec![], Some((2.min(chain_len), "=O".to_string())))),
        "anoic acid" => Some((vec![], Some((1, "=O)(O".to_string())))),
        "anamide" => Some((vec![], Some((1, "=O)(N".to_string())))),
        "anamine" => Some((vec![], Some((1, "N".to_string())))),
        "anethiol" => Some((vec![], Some((1, "S".to_string())))),
        "anenitrile" | "anonitrile" => Some((vec![], Some((1, "#N".to_string())))),
        _ => {
            // Handle infix locants like "an-2-ol", "an-3-one", "-2-ene", "-1-yne", "an-2-amine"
            let clean = s.strip_prefix("an").unwrap_or(s);
            let parts: Vec<&str> = clean.split('-').filter(|p| !p.is_empty()).collect();
            if parts.len() == 2 {
                if let Ok(loc) = parts[0].parse::<usize>() {
                    let pos = loc.saturating_sub(1);
                    return match parts[1] {
                        "ene" => Some((vec![(pos, 2)], None)),
                        "yne" => Some((vec![(pos, 3)], None)),
                        "ol" => Some((vec![], Some((loc, "O".to_string())))),
                        "one" => Some((vec![], Some((loc, "=O".to_string())))),
                        "amine" => Some((vec![], Some((loc, "N".to_string())))),
                        "thiol" => Some((vec![], Some((loc, "S".to_string())))),
                        _ => None,
                    };
                }
            }
            None
        }
    }
}

/// Parses substituent prefixes like `"2,4-dichloro-3-methyl-"` into `Some(Vec<(locant, smiles)>)`,
/// returning `None` if a non-empty prefix contains unrecognized substituent tokens.
fn parse_substituent_prefix_checked(prefix: &str) -> Option<Vec<(usize, String)>> {
    let mut out = Vec::new();
    let norm = prefix
        .trim()
        .trim_end_matches('-')
        .replace("tert-butyl", "tertbutyl")
        .replace("t-butyl", "tertbutyl");
    if norm.is_empty() {
        return Some(out);
    }
    let tokens: Vec<&str> = norm.split('-').filter(|t| !t.is_empty()).collect();
    let mut pending_locs: Vec<usize> = Vec::new();
    for tok in tokens {
        let nums: Vec<usize> = tok.split(',').filter_map(|x| x.parse::<usize>().ok()).collect();
        if !nums.is_empty() && nums.len() == tok.split(',').count() {
            pending_locs = nums;
            continue;
        }
        let mut rem = tok;
        while !rem.is_empty() {
            let stripped = strip_multiplier(rem);
            if let Some(smi) = subst_to_smiles(stripped) {
                let locs = if pending_locs.is_empty() { vec![1] } else { std::mem::take(&mut pending_locs) };
                for l in locs {
                    out.push((l, smi.to_string()));
                }
                break;
            }
            let mut found_prefix = None;
            for cand in [
                "trifluoromethyl", "tertbutyl", "isopropyl", "methoxy", "ethoxy", "hydroxy",
                "fluoro", "chloro", "bromo", "iodo", "amino", "nitro", "cyano", "methyl", "ethyl",
                "propyl", "butyl", "phenyl", "acetyl", "formyl", "carboxy", "oxo",
            ] {
                let base = strip_multiplier(rem);
                if let Some(rest) = base.strip_prefix(cand) {
                    found_prefix = Some((cand, rest));
                    break;
                }
            }
            if let Some((cand, rest)) = found_prefix {
                if let Some(smi) = subst_to_smiles(cand) {
                    let locs = if pending_locs.is_empty() { vec![1] } else { std::mem::take(&mut pending_locs) };
                    for l in locs {
                        out.push((l, smi.to_string()));
                    }
                    rem = rest;
                } else {
                    return None;
                }
            } else {
                return None;
            }
        }
    }
    if !pending_locs.is_empty() {
        return None;
    }
    Some(out)
}

/// Returns `true` when `molecule_to_name` could only produce a generic formula fallback.
pub fn is_fallback_chemical_name(name: &str) -> bool {
    let t = name.trim();
    t.is_empty() || t.starts_with("Compound (")
}

/// Converts a `ChemMolecule` into its IUPAC or canonical preferred chemical name!
pub fn molecule_to_name(mol: &ChemMolecule, smiles: &str, dict: &NameDictionary) -> String {
    if mol.atoms.is_empty() {
        return String::new();
    }
    // 1. Check user dictionary by exact SMILES match
    for (name, smi) in &dict.user_entries {
        if smi == smiles {
            return name.clone();
        }
    }
    let norm_smi = smiles.trim();
    let formula = mol.compute_formula_properties().get_formula();
    let arom = CcAromaticityInfo::perceive(mol);
    let val = CcVal::compute(mol);
    let heavy_n = mol.atoms.iter().filter(|a| a.z != 1).count();
    let arom_c: Vec<usize> = (0..mol.atoms.len())
        .filter(|&i| arom.is_atom_aromatic(i) && mol.atoms[i].z == 6)
        .collect();
    let arom_n: Vec<usize> = (0..mol.atoms.len())
        .filter(|&i| arom.is_atom_aromatic(i) && mol.atoms[i].z == 7)
        .collect();
    let arom_o: Vec<usize> = (0..mol.atoms.len())
        .filter(|&i| arom.is_atom_aromatic(i) && mol.atoms[i].z == 8)
        .collect();
    let arom_s: Vec<usize> = (0..mol.atoms.len())
        .filter(|&i| arom.is_atom_aromatic(i) && mol.atoms[i].z == 16)
        .collect();

    // Fast exact recognition of flagship drugs, solvents, bridged cages, amino acids & parent rings
    let c_deg3 = (0..mol.atoms.len())
        .filter(|&i| mol.atoms[i].z == 6 && mol.neighbors(i).iter().filter(|&&(nb, _)| mol.atoms[nb].z != 1).count() == 3)
        .count();
    let c_deg2 = (0..mol.atoms.len())
        .filter(|&i| mol.atoms[i].z == 6 && mol.neighbors(i).iter().filter(|&&(nb, _)| mol.atoms[nb].z != 1).count() == 2)
        .count();
    let heavy_bonds = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom < mol.atoms.len() && b.end_atom < mol.atoms.len() && mol.atoms[b.begin_atom].z != 1 && mol.atoms[b.end_atom].z != 1)
        .count();
    if arom_c.is_empty() && arom_n.is_empty() && mol.bonds.iter().all(|b| b.order == 1) {
        if formula == "C10H16" && heavy_n == 10 && heavy_bonds == 12 && c_deg3 == 4 && c_deg2 == 6 {
            return "adamantane".to_string();
        }
        if formula == "C7H12" && heavy_n == 7 && heavy_bonds == 8 && c_deg3 == 2 && c_deg2 == 5 {
            return "norbornane (bicyclo[2.2.1]heptane)".to_string();
        }
        if formula == "C8H8" && heavy_n == 8 && heavy_bonds == 12 && c_deg3 == 8 {
            return "cubane".to_string();
        }
        if formula == "C8H14" && heavy_n == 8 && heavy_bonds == 9 && c_deg3 == 2 && c_deg2 == 6 {
            return "bicyclo[2.2.2]octane".to_string();
        }
        if formula == "C10H18" && heavy_n == 10 && heavy_bonds == 11 && c_deg3 == 2 && c_deg2 == 8 {
            return "decahydronaphthalene (decalin)".to_string();
        }
        if formula == "C12H24O6" && heavy_n == 18 && heavy_bonds == 18 {
            return "18-crown-6".to_string();
        }
    }
    if formula == "C20H14N4" && heavy_n == 24 && heavy_bonds == 28 {
        return "porphine".to_string();
    }

    match (formula.as_str(), arom_c.len(), arom_n.len(), heavy_n) {
        ("C9H8O4", 6, 0, 13) => return "2-acetoxybenzoic acid (aspirin)".to_string(),
        ("C8H10N4O2", 5, 4, 14) => return "1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione (caffeine)".to_string(),
        ("C13H18O2", 6, 0, 15) => return "2-(4-isobutylphenyl)propanoic acid (ibuprofen)".to_string(),
        ("C8H9NO2", 6, 0, 11) => return "N-(4-hydroxyphenyl)acetamide (acetaminophen)".to_string(),
        ("C14H14O3", 10, 0, 17) => return "2-(6-methoxynaphthalen-2-yl)propanoic acid (naproxen)".to_string(),
        ("C17H19NO3", 6, 0, 21) => return "morphine".to_string(),
        ("C10H14N2", 5, 1, 12) => return "3-(1-methylpyrrolidin-2-yl)pyridine (nicotine)".to_string(),
        ("C8H11NO2", 6, 0, 11) => return "4-(2-aminoethyl)benzene-1,2-diol (dopamine)".to_string(),
        ("C10H12N2O", 8, 1, 13) => return "3-(2-aminoethyl)-1H-indol-5-ol (serotonin)".to_string(),
        ("C10H8", 10, 0, 10) => return "naphthalene".to_string(),
        ("C14H10", 14, 0, 14) => return "anthracene".to_string(),
        ("C12H10", 12, 0, 12) => return "1,1'-biphenyl".to_string(),
        ("C8H7N", 8, 1, 9) => return "1H-indole".to_string(),
        ("C9H7N", 9, 1, 10) => return "quinoline".to_string(),
        ("C5H4N4", 5, 4, 9) => return "9H-purine".to_string(),
        ("C6H6", 6, 0, 6) => return "benzene".to_string(),
        ("C5H5N", 5, 1, 6) => return "pyridine".to_string(),
        ("C4H4N2", 4, 2, 6) => return "pyrimidine".to_string(),
        ("C4H5N", 4, 1, 5) => return "1H-pyrrole".to_string(),
        ("C3H4N2", 3, 2, 5) => return "1H-imidazole".to_string(),
        ("C4H4O", 4, 0, 5) => return "furan".to_string(),
        ("C4H4S", 4, 0, 5) => return "thiophene".to_string(),
        ("C3H3NO", 3, 1, 5) => return "1,3-oxazole".to_string(),
        ("C3H3NS", 3, 1, 5) => return "1,3-thiazole".to_string(),
        ("C7H8", 6, 0, 7) => return "toluene".to_string(),
        ("C6H6O", 6, 0, 7) => return "phenol".to_string(),
        ("C6H7N", 6, 0, 7) => return "aniline".to_string(),
        ("C7H8O", 6, 0, 8) if mol.atoms.iter().any(|a| a.z == 8 && mol.neighbors(mol.atoms.iter().position(|x| x.id == a.id).unwrap_or(0)).len() >= 2) => {
            return "methoxybenzene (anisole)".to_string();
        }
        ("C7H6O2", 6, 0, 9) => return "benzoic acid".to_string(),
        ("C7H6O", 6, 0, 8) => return "benzaldehyde".to_string(),
        ("C8H8", 6, 0, 8) => return "styrene".to_string(),
        ("C7H5N", 6, 0, 8) => return "benzonitrile".to_string(),
        ("C6H5NO2", 6, 0, 9) => return "nitrobenzene".to_string(),
        ("C5H11N", 0, 0, 6) if mol.bonds.len() == 6 => return "piperidine".to_string(),
        ("C4H10N2", 0, 0, 6) if mol.bonds.len() == 6 => return "piperazine".to_string(),
        ("C4H9NO", 0, 0, 6) if mol.bonds.len() == 6 => return "morpholine".to_string(),
        ("C4H9N", 0, 0, 5) if mol.bonds.len() == 5 => return "pyrrolidine".to_string(),
        ("C4H8O", 0, 0, 5) if mol.bonds.len() == 5 => return "oxolane (THF)".to_string(),
        ("C4H8O2", 0, 0, 6) if mol.bonds.len() == 6 => return "1,4-dioxane".to_string(),
        ("H2O", 0, 0, 1) => return "water".to_string(),
        ("H3N", 0, 0, 1) => return "ammonia".to_string(),
        ("CH4", 0, 0, 1) => return "methane".to_string(),
        ("CH4O", 0, 0, 2) => return "methanol".to_string(),
        ("CH2O", 0, 0, 2) => return "methanal (formaldehyde)".to_string(),
        ("CH2O2", 0, 0, 3) => return "formic acid".to_string(),
        ("CHCl3", 0, 0, 4) => return "trichloromethane (chloroform)".to_string(),
        ("CH2Cl2", 0, 0, 3) => return "dichloromethane".to_string(),
        ("CH4N2O", 0, 0, 4) => return "urea".to_string(),
        ("C2H6OS", 0, 0, 4) => return "methylsulfinylmethane (DMSO)".to_string(),
        ("C3H7NO", 0, 0, 5) if norm_smi.contains("CN(C)C=O") => return "N,N-dimethylformamide (DMF)".to_string(),
        ("C2H6O", 0, 0, 3) if norm_smi.contains("CCO") || norm_smi.contains("OCC") => return "ethanol".to_string(),
        ("C2H6O", 0, 0, 3) if norm_smi == "COC" => return "methoxymethane (dimethyl ether)".to_string(),
        ("C2H4O", 0, 0, 3) if mol.bonds.iter().any(|b| b.order == 2) => return "ethanal (acetaldehyde)".to_string(),
        ("C2H3N", 0, 0, 3) if mol.bonds.iter().any(|b| b.order == 3) => return "acetonitrile".to_string(),
        ("C2H4O2", 0, 0, 4) if mol.atoms.iter().filter(|a| a.z == 8).count() == 2 && mol.bonds.iter().any(|b| b.order == 2) => {
            if let Some(est) = name_simple_ester(mol) {
                return est;
            }
            return "acetic acid".to_string();
        }
        ("C3H6O", 0, 0, 4) if mol.bonds.iter().any(|b| b.order == 2) && mol.atoms.iter().any(|a| a.z == 8) => {
            let is_ketone = mol.atoms.iter().enumerate().any(|(i, a)| {
                a.z == 6
                    && mol.neighbors(i).iter().any(|&(nb, bi)| mol.atoms[nb].z == 8 && mol.bonds[bi].order == 2)
                    && mol.neighbors(i).iter().filter(|&&(nb, _)| mol.atoms[nb].z == 6).count() == 2
            });
            return if is_ketone { "propan-2-one (acetone)".to_string() } else { "propanal".to_string() };
        }
        _ => {}
    }

    // 2. Simple esters (`alkyl alkanoate` / `alkyl benzoate`, e.g. ethyl acetate, methyl benzoate)
    if let Some(est_name) = name_simple_ester(mol) {
        return prepend_stereo(mol, &est_name);
    }

    // 3. Substituted benzene ring systematic naming
    if arom_c.len() == 6 && arom_n.is_empty() && arom_o.is_empty() && arom_s.is_empty() {
        if let Some(name) = name_substituted_benzene(mol, &arom, &arom_c) {
            return prepend_stereo(mol, &name);
        }
    }

    // 4. Substituted aromatic heterocycles (pyridine, furan, thiophene, pyrrole)
    if let Some(het_name) = name_substituted_heterocycle(mol, &arom, &arom_c, &arom_n, &arom_o, &arom_s) {
        return prepend_stereo(mol, &het_name);
    }

    // 5. Substituted & simple monocyclic cycloalkanes / cycloalkenes / cycloalkanols / cycloalkanones
    let carbon_indices: Vec<usize> = (0..mol.atoms.len()).filter(|&i| mol.atoms[i].z == 6).collect();
    if let Some(ring_name) = name_substituted_cycloalkane(mol) {
        return prepend_stereo(mol, &ring_name);
    }

    // 6. Simple acyclic ethers (`R-O-R'`, e.g. diethyl ether, 1-methoxypropane)
    if let Some(ether_name) = name_simple_ether(mol) {
        return prepend_stereo(mol, &ether_name);
    }

    // 7. Systematic acyclic carbon-chain IUPAC naming (longest chain + principal group + lowest-locant substituents)
    if let Some(chain_name) = name_acyclic_chain(mol, &val, &carbon_indices) {
        return prepend_stereo(mol, &chain_name);
    }

    // Fallback: Hill formula descriptor (triggers optional live PubChem lookup!)
    format!("Compound ({formula})")
}

fn prepend_stereo(mol: &ChemMolecule, base: &str) -> String {
    let st = perceive_stereochemistry(mol);
    let mut tags: Vec<(usize, &'static str)> = Vec::new();
    for (&ai, &code) in &st.atom_stereo {
        match code {
            2 => tags.push((ai, "R")),
            3 => tags.push((ai, "S")),
            4 => tags.push((ai, "r")),
            5 => tags.push((ai, "s")),
            _ => {}
        }
    }
    for (&bi, &code) in &st.bond_stereo {
        match code {
            10 => tags.push((1000 + bi, "Z")),
            11 => tags.push((1000 + bi, "E")),
            _ => {}
        }
    }
    tags.sort_by_key(|&(k, _)| k);
    if tags.is_empty() {
        base.to_string()
    } else if tags.len() == 1 {
        format!("({})-{base}", tags[0].1)
    } else {
        let joined: Vec<&str> = tags.iter().map(|&(_, t)| t).collect();
        format!("({})-{base}", joined.join(","))
    }
}

/// Names simple carboxylic esters `R1-C(=O)-O-R2` (alkyl alkanoates & alkyl benzoates).
fn name_simple_ester(mol: &ChemMolecule) -> Option<String> {
    // Find bridging oxygen bonded to exactly 2 carbons (no H): one carbonyl C(=O) and one non-carbonyl C
    for (oi, a) in mol.atoms.iter().enumerate() {
        if a.z != 8 {
            continue;
        }
        let nbrs: Vec<(usize, usize)> = mol
            .neighbors(oi)
            .into_iter()
            .filter(|&(nb, _)| mol.atoms[nb].z != 1)
            .collect();
        if nbrs.len() != 2 || nbrs.iter().any(|&(nb, bi)| mol.atoms[nb].z != 6 || mol.bonds[bi].order != 1) {
            continue;
        }
        let c0 = nbrs[0].0;
        let c1 = nbrs[1].0;
        let is_carbonyl = |ci: usize| {
            mol.neighbors(ci)
                .iter()
                .any(|&(nb, bi)| nb != oi && mol.atoms[nb].z == 8 && mol.bonds[bi].order == 2)
        };
        let (acyl_c, alk_c) = if is_carbonyl(c0) && !is_carbonyl(c1) {
            (c0, c1)
        } else if is_carbonyl(c1) && !is_carbonyl(c0) {
            (c1, c0)
        } else {
            continue;
        };
        // Ensure bridging bond splits molecule cleanly (ester is acyclic bridge, not lactone)
        let acyl_sub = collect_subtree(mol, acyl_c, oi);
        let alk_sub = collect_subtree(mol, alk_c, oi);
        if acyl_sub.iter().any(|x| alk_sub.contains(x)) {
            continue; // cyclic lactone
        }
        // Only handle pure C/H/O simple esters here (leave complex drug esters like aspirin to higher rules/PubChem)
        if acyl_sub.iter().chain(alk_sub.iter()).any(|&i| !matches!(mol.atoms[i].z, 1 | 6 | 8)) {
            continue;
        }
        let alk_c_cnt = alk_sub.iter().filter(|&&i| mol.atoms[i].z == 6).count();
        let alk_o_cnt = alk_sub.iter().filter(|&&i| mol.atoms[i].z == 8).count();
        if alk_o_cnt > 0 || alk_c_cnt == 0 {
            continue;
        }
        let arom = CcAromaticityInfo::perceive(mol);
        let alk_name = if alk_c_cnt == 6 && alk_sub.iter().all(|&i| mol.atoms[i].z != 6 || arom.is_atom_aromatic(i)) {
            "phenyl".to_string()
        } else if let Some(&(r, _)) = CHAIN_ROOTS.iter().find(|&&(_, n)| n == alk_c_cnt) {
            format!("{r}yl")
        } else {
            continue;
        };
        let acyl_c_cnt = acyl_sub.iter().filter(|&&i| mol.atoms[i].z == 6).count();
        let acyl_o_cnt = acyl_sub.iter().filter(|&&i| mol.atoms[i].z == 8).count();
        if acyl_o_cnt != 1 {
            continue; // only the single =O on acyl side
        }
        let acid_part = if acyl_c_cnt == 7 && acyl_sub.iter().filter(|&&i| arom.is_atom_aromatic(i)).count() == 6 {
            "benzoate".to_string()
        } else if acyl_c_cnt == 1 {
            "formate".to_string()
        } else if acyl_c_cnt == 2 {
            "acetate".to_string()
        } else if let Some(&(r, _)) = CHAIN_ROOTS.iter().find(|&&(_, n)| n == acyl_c_cnt) {
            format!("{r}anoate")
        } else {
            continue;
        };
        return Some(format!("{alk_name} {acid_part}"));
    }
    None
}

/// Names simple acyclic ethers `R1-O-R2` (when molecule contains exactly one oxygen, no rings/carbonyls).
fn name_simple_ether(mol: &ChemMolecule) -> Option<String> {
    let heavy: Vec<usize> = (0..mol.atoms.len()).filter(|&i| mol.atoms[i].z != 1).collect();
    if heavy.iter().any(|&i| !matches!(mol.atoms[i].z, 6 | 8)) {
        return None;
    }
    let oxygens: Vec<usize> = heavy.iter().copied().filter(|&i| mol.atoms[i].z == 8).collect();
    if oxygens.len() != 1 || mol.bonds.iter().any(|b| b.order != 1) {
        return None;
    }
    let oi = oxygens[0];
    let nbrs: Vec<usize> = mol
        .neighbors(oi)
        .into_iter()
        .map(|(nb, _)| nb)
        .filter(|&nb| mol.atoms[nb].z == 6)
        .collect();
    if nbrs.len() != 2 {
        return None;
    }
    let s1 = collect_subtree(mol, nbrs[0], oi);
    let s2 = collect_subtree(mol, nbrs[1], oi);
    if s1.iter().any(|x| s2.contains(x)) {
        return None; // cyclic ether
    }
    let n1 = s1.iter().filter(|&&i| mol.atoms[i].z == 6).count();
    let n2 = s2.iter().filter(|&&i| mol.atoms[i].z == 6).count();
    let (small, large) = if n1 <= n2 { (n1, n2) } else { (n2, n1) };
    let r_small = CHAIN_ROOTS.iter().find(|&&(_, n)| n == small)?.0;
    let r_large = CHAIN_ROOTS.iter().find(|&&(_, n)| n == large)?.0;
    if small == large {
        Some(format!("{r_small}oxy{r_large}ane (di{r_small}yl ether)"))
    } else {
        Some(format!("{r_small}oxy{r_large}ane"))
    }
}

fn collect_subtree(mol: &ChemMolecule, start: usize, blocked: usize) -> Vec<usize> {
    let mut visited = vec![false; mol.atoms.len()];
    visited[blocked] = true;
    visited[start] = true;
    let mut stack = vec![start];
    let mut out = Vec::new();
    while let Some(u) = stack.pop() {
        out.push(u);
        for &(v, _) in &mol.neighbors(u) {
            if !visited[v] && mol.atoms[v].z != 1 {
                visited[v] = true;
                stack.push(v);
            }
        }
    }
    out
}

/// Names substituted 5- and 6-membered aromatic heterocycles (`pyridine`, `furan`, `thiophene`, `1H-pyrrole`).
fn name_substituted_heterocycle(
    mol: &ChemMolecule,
    arom: &CcAromaticityInfo,
    arom_c: &[usize],
    arom_n: &[usize],
    arom_o: &[usize],
    arom_s: &[usize],
) -> Option<String> {
    let (hetero_idx, ring_len, parent) = if arom_c.len() == 5 && arom_n.len() == 1 && arom_o.is_empty() && arom_s.is_empty() {
        (arom_n[0], 6usize, "pyridine")
    } else if arom_c.len() == 4 && arom_o.len() == 1 && arom_n.is_empty() && arom_s.is_empty() {
        (arom_o[0], 5usize, "furan")
    } else if arom_c.len() == 4 && arom_s.len() == 1 && arom_n.is_empty() && arom_o.is_empty() {
        (arom_s[0], 5usize, "thiophene")
    } else if arom_c.len() == 4 && arom_n.len() == 1 && arom_o.is_empty() && arom_s.is_empty() {
        (arom_n[0], 5usize, "1H-pyrrole")
    } else {
        return None;
    };
    let mut ring = vec![hetero_idx];
    while ring.len() < ring_len {
        let cur = *ring.last().unwrap();
        let next = mol
            .neighbors(cur)
            .into_iter()
            .map(|(nb, _)| nb)
            .find(|&nb| arom.is_atom_aromatic(nb) && !ring.contains(&nb))?;
        ring.push(next);
    }
    let mut pos_subs: Vec<(usize, String)> = Vec::new();
    for (idx, &ri) in ring.iter().enumerate() {
        for &(nb, _) in &mol.neighbors(ri) {
            if !arom.is_atom_aromatic(nb) && mol.atoms[nb].z != 1 {
                pos_subs.push((idx, classify_substituent(mol, nb, ri)));
            }
        }
    }
    if pos_subs.is_empty() {
        return Some(parent.to_string());
    }
    // Heteroatom is always locant 1; choose CW vs CCW for lowest substituent locants
    let mut cw: Vec<(usize, String)> = pos_subs.iter().map(|(i, s)| (i + 1, s.clone())).collect();
    let mut ccw: Vec<(usize, String)> = pos_subs
        .iter()
        .map(|(i, s)| (if *i == 0 { 1 } else { ring_len + 1 - i }, s.clone()))
        .collect();
    let mut cw_locs: Vec<usize> = cw.iter().map(|x| x.0).collect();
    let mut ccw_locs: Vec<usize> = ccw.iter().map(|x| x.0).collect();
    cw_locs.sort_unstable();
    ccw_locs.sort_unstable();
    let chosen = if ccw_locs < cw_locs { &mut ccw } else { &mut cw };
    Some(format!("{}{}", format_substituents(chosen), parent))
}

fn classify_substituent(mol: &ChemMolecule, sub_root: usize, parent_atom: usize) -> String {
    let z = mol.atoms[sub_root].z;
    let nbrs: Vec<(usize, usize)> = mol
        .neighbors(sub_root)
        .into_iter()
        .filter(|&(nb, _)| nb != parent_atom && mol.atoms[nb].z != 1)
        .collect();
    match z {
        9 => "fluoro".to_string(),
        17 => "chloro".to_string(),
        35 => "bromo".to_string(),
        53 => "iodo".to_string(),
        8 => {
            if nbrs.is_empty() {
                "hydroxy".to_string()
            } else if nbrs.len() == 1 && mol.atoms[nbrs[0].0].z == 6 {
                let c_sub = count_subtree_carbons(mol, nbrs[0].0, sub_root);
                match c_sub {
                    1 => "methoxy".to_string(),
                    2 => "ethoxy".to_string(),
                    _ => "alkoxy".to_string(),
                }
            } else {
                "oxy".to_string()
            }
        }
        7 => {
            let o_cnt = nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z == 8).count();
            if o_cnt >= 2 {
                "nitro".to_string()
            } else if let Some(&(c_idx, _)) = nbrs.iter().find(|&&(nb, _)| {
                mol.atoms[nb].z == 6
                    && mol.neighbors(nb).iter().any(|&(k, bi)| mol.atoms[k].z == 8 && mol.bonds[bi].order == 2)
            }) {
                let c_sub = count_subtree_carbons(mol, c_idx, sub_root);
                if c_sub == 2 {
                    "acetamido".to_string()
                } else {
                    "amido".to_string()
                }
            } else {
                "amino".to_string()
            }
        }
        16 => "sulfanyl".to_string(),
        6 => {
            let has_dbl_o = nbrs.iter().any(|&(nb, bi)| mol.atoms[nb].z == 8 && mol.bonds[bi].order == 2);
            let has_sgl_o = nbrs.iter().any(|&(nb, bi)| mol.atoms[nb].z == 8 && mol.bonds[bi].order == 1);
            let has_trp_n = nbrs.iter().any(|&(nb, bi)| mol.atoms[nb].z == 7 && mol.bonds[bi].order == 3);
            if has_dbl_o && has_sgl_o {
                "carboxy".to_string()
            } else if has_dbl_o && nbrs.len() == 1 {
                "formyl".to_string()
            } else if has_trp_n {
                "cyano".to_string()
            } else {
                let cnt = count_subtree_carbons(mol, sub_root, parent_atom);
                match cnt {
                    1 => "methyl".to_string(),
                    2 => "ethyl".to_string(),
                    3 => "propyl".to_string(),
                    4 => "butyl".to_string(),
                    _ => "alkyl".to_string(),
                }
            }
        }
        _ => "substituted".to_string(),
    }
}

fn count_subtree_carbons(mol: &ChemMolecule, start: usize, blocked: usize) -> usize {
    let mut visited = vec![false; mol.atoms.len()];
    visited[blocked] = true;
    visited[start] = true;
    let mut stack = vec![start];
    let mut c_cnt = 0usize;
    while let Some(u) = stack.pop() {
        if mol.atoms[u].z == 6 {
            c_cnt += 1;
        }
        for &(v, _) in &mol.neighbors(u) {
            if !visited[v] && mol.atoms[v].z != 1 {
                visited[v] = true;
                stack.push(v);
            }
        }
    }
    c_cnt
}

fn name_substituted_benzene(
    mol: &ChemMolecule,
    arom: &CcAromaticityInfo,
    arom_c: &[usize],
) -> Option<String> {
    // Walk ring in order
    let mut ring = vec![arom_c[0]];
    while ring.len() < 6 {
        let cur = *ring.last().unwrap();
        let next = mol
            .neighbors(cur)
            .into_iter()
            .map(|(nb, _)| nb)
            .find(|&nb| arom.is_atom_aromatic(nb) && mol.atoms[nb].z == 6 && !ring.contains(&nb))?;
        ring.push(next);
    }
    // Collect substituents at each ring position 0..6
    let mut pos_subs: Vec<(usize, String)> = Vec::new();
    for (idx, &ci) in ring.iter().enumerate() {
        for &(nb, _) in &mol.neighbors(ci) {
            if !arom.is_atom_aromatic(nb) && mol.atoms[nb].z != 1 {
                pos_subs.push((idx, classify_substituent(mol, nb, ci)));
            }
        }
    }
    if pos_subs.is_empty() {
        return Some("benzene".to_string());
    }
    // Check principal parent group priority: carboxy -> benzoic acid, formyl -> benzaldehyde, hydroxy -> phenol, amino -> aniline, methyl -> toluene
    let parent_priority = [
        ("carboxy", "benzoic acid"),
        ("formyl", "benzaldehyde"),
        ("cyano", "benzonitrile"),
        ("acetamido", "acetanilide"),
        ("hydroxy", "phenol"),
        ("amino", "aniline"),
        ("methyl", "toluene"),
    ];
    let mut parent_name = "benzene";
    let mut anchor_idx = pos_subs[0].0;
    let mut consumed_idx: Option<usize> = None;
    for &(sub_key, pname) in &parent_priority {
        if let Some((si, (ridx, _))) = pos_subs.iter().enumerate().find(|(_, (_, s))| s == sub_key) {
            parent_name = pname;
            anchor_idx = *ridx;
            consumed_idx = Some(si);
            break;
        }
    }
    let raw_rem: Vec<(usize, String)> = pos_subs
        .into_iter()
        .enumerate()
        .filter(|&(si, _)| Some(si) != consumed_idx)
        .map(|(_, (ridx, s))| (ridx, s))
        .collect();
    if raw_rem.is_empty() {
        return Some(parent_name.to_string());
    }
    // Choose a SINGLE consistent CW or CCW ring numbering (and best anchor for bare benzene)
    let anchors: Vec<usize> = if parent_name == "benzene" {
        raw_rem.iter().map(|(r, _)| *r).collect()
    } else {
        vec![anchor_idx]
    };
    let mut best_locs: Option<Vec<usize>> = None;
    let mut best_rem: Vec<(usize, String)> = Vec::new();
    for anc in anchors {
        for dir in [1isize, -1isize] {
            let mapped: Vec<(usize, String)> = raw_rem
                .iter()
                .map(|(ridx, s)| {
                    let loc = if dir == 1 {
                        (*ridx + 6 - anc) % 6 + 1
                    } else {
                        (anc + 6 - *ridx) % 6 + 1
                    };
                    (loc, s.clone())
                })
                .collect();
            let mut locs: Vec<usize> = mapped.iter().map(|x| x.0).collect();
            locs.sort_unstable();
            if best_locs.as_ref().map(|b| locs < *b).unwrap_or(true) {
                best_locs = Some(locs);
                best_rem = mapped;
            }
        }
    }
    if best_rem.len() == 1 && parent_name == "benzene" {
        return Some(format!("{}{}", best_rem[0].1, parent_name));
    }
    Some(format!("{}{}", format_substituents(&best_rem), parent_name))
}

/// Names monocyclic 3..=8-membered cycloalkanes, cycloalkenes, cycloalkanols, cycloalkanones,
/// and their substituted derivatives using IUPAC lowest-locant ring numbering!
fn name_substituted_cycloalkane(mol: &ChemMolecule) -> Option<String> {
    // Find a single non-aromatic carbon ring (size 3..=8) via DFS cycle search
    let n_atoms = mol.atoms.len();
    let mut found_ring: Option<Vec<usize>> = None;
    for start in 0..n_atoms {
        if mol.atoms[start].z != 6 {
            continue;
        }
        let mut stack: Vec<(usize, Vec<usize>)> = vec![(start, vec![start])];
        while let Some((u, path)) = stack.pop() {
            for &(v, _) in &mol.neighbors(u) {
                if mol.atoms[v].z != 6 {
                    continue;
                }
                if v == start && path.len() >= 3 && path.len() <= 8 {
                    if found_ring.as_ref().map(|r| path.len() > r.len()).unwrap_or(true) {
                        found_ring = Some(path.clone());
                    }
                } else if !path.contains(&v) && path.len() < 8 && v > start {
                    let mut next = path.clone();
                    next.push(v);
                    stack.push((v, next));
                }
            }
        }
    }
    let ring = found_ring?;
    let rlen = ring.len();
    let in_ring: HashMap<usize, usize> = ring.iter().enumerate().map(|(i, &a)| (a, i)).collect();

    // Ensure there are no OTHER cycles (fused/spiro systems fall through to PubChem/flagship)
    let heavy_n = (0..n_atoms).filter(|&i| mol.atoms[i].z != 1).count();
    let heavy_b = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom < n_atoms && b.end_atom < n_atoms && mol.atoms[b.begin_atom].z != 1 && mol.atoms[b.end_atom].z != 1)
        .count();
    if heavy_b != heavy_n {
        return None; // monocyclic graph always satisfies E == V
    }

    let root = CHAIN_ROOTS.iter().find(|&&(_, len)| len == rlen)?.0;
    // Detect ring C=C double bonds & exocyclic groups
    let mut ring_dbl: Option<usize> = None;
    for i in 0..rlen {
        let a = ring[i];
        let b_atom = ring[(i + 1) % rlen];
        if mol.bonds.iter().any(|bd| {
            ((bd.begin_atom == a && bd.end_atom == b_atom) || (bd.begin_atom == b_atom && bd.end_atom == a))
                && bd.order == 2
        }) {
            ring_dbl = Some(i);
            break;
        }
    }

    let mut pos_subs: Vec<(usize, String, bool)> = Vec::new(); // (ring_idx, label, is_ketone)
    for (idx, &ci) in ring.iter().enumerate() {
        for &(nb, bi) in &mol.neighbors(ci) {
            if in_ring.contains_key(&nb) || mol.atoms[nb].z == 1 {
                continue;
            }
            if mol.atoms[nb].z == 8 && mol.bonds[bi].order == 2 {
                pos_subs.push((idx, "one".to_string(), true));
            } else {
                pos_subs.push((idx, classify_substituent(mol, nb, ci), false));
            }
        }
    }

    // Principal group on ring: ketone ("anone") > alcohol ("anol") > amine ("anamine") > alkene ("ene") > alkane ("ane")
    let (suffix, anchor_opt, consumed_idx) = if let Some((si, (ridx, _, _))) =
        pos_subs.iter().enumerate().find(|(_, (_, _, is_k))| *is_k)
    {
        ("anone", Some(*ridx), Some(si))
    } else if let Some((si, (ridx, _, _))) =
        pos_subs.iter().enumerate().find(|(_, (_, s, _))| s == "hydroxy")
    {
        ("anol", Some(*ridx), Some(si))
    } else if let Some((si, (ridx, _, _))) =
        pos_subs.iter().enumerate().find(|(_, (_, s, _))| s == "amino")
    {
        ("anamine", Some(*ridx), Some(si))
    } else if let Some(d_idx) = ring_dbl {
        ("ene", Some(d_idx), None)
    } else {
        ("ane", None, None)
    };

    let raw_rem: Vec<(usize, String)> = pos_subs
        .into_iter()
        .enumerate()
        .filter(|&(si, _)| Some(si) != consumed_idx)
        .map(|(_, (ridx, s, _))| (ridx, s))
        .collect();

    if raw_rem.is_empty() {
        return Some(format!("cyclo{root}{suffix}"));
    }

    // Pick lowest-locant ring orientation (for cycloalkenes, C1=C2 must traverse the double bond!)
    let cand_pairs: Vec<(usize, isize)> = if suffix == "ene" {
        if let Some(d_idx) = ring_dbl {
            vec![(d_idx, 1), ((d_idx + 1) % rlen, -1)]
        } else {
            vec![(0, 1), (0, -1)]
        }
    } else if let Some(a) = anchor_opt {
        vec![(a, 1), (a, -1)]
    } else {
        raw_rem.iter().flat_map(|&(r, _)| [(r, 1), (r, -1)]).collect()
    };
    let mut best_locs: Option<Vec<usize>> = None;
    let mut best_rem: Vec<(usize, String)> = Vec::new();
    for (anc, dir) in cand_pairs {
        let mapped: Vec<(usize, String)> = raw_rem
            .iter()
            .map(|(r, s)| {
                let loc = if dir == 1 {
                    (*r + rlen - anc) % rlen + 1
                } else {
                    (anc + rlen - *r) % rlen + 1
                };
                (loc, s.clone())
            })
            .collect();
        let mut locs: Vec<usize> = mapped.iter().map(|x| x.0).collect();
        locs.sort_unstable();
        if best_locs.as_ref().map(|b| locs < *b).unwrap_or(true) {
            best_locs = Some(locs);
            best_rem = mapped;
        }
    }
    if best_rem.len() == 1 && suffix == "ane" {
        return Some(format!("{}cyclo{root}{suffix}", best_rem[0].1));
    }
    Some(format!("{}cyclo{root}{suffix}", format_substituents(&best_rem)))
}

fn name_acyclic_chain(mol: &ChemMolecule, _val: &CcVal, carbons: &[usize]) -> Option<String> {
    if carbons.is_empty() {
        return None;
    }
    // Reject cyclic graphs up front (avoids exponential simple-path search on polycycles/cages)
    let n_atoms = mol.atoms.len();
    let heavy_n = (0..n_atoms).filter(|&i| mol.atoms[i].z != 1).count();
    let heavy_b = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom < n_atoms && b.end_atom < n_atoms && b.begin_atom != b.end_atom && mol.atoms[b.begin_atom].z != 1 && mol.atoms[b.end_atom].z != 1)
        .count();
    if heavy_b >= heavy_n {
        return None;
    }
    // Find longest simple carbon chain via bounded DFS
    let mut best_chain: Vec<usize> = Vec::new();
    let mut budget = 10_000usize;
    for &start in carbons {
        if budget == 0 {
            break;
        }
        let mut path = vec![start];
        let mut visited = vec![false; n_atoms];
        visited[start] = true;
        dfs_longest_c_chain(mol, start, &mut visited, &mut path, &mut best_chain, &mut budget);
    }
    let n = best_chain.len();
    if n == 0 || n > 12 {
        return None;
    }
    let root = CHAIN_ROOTS.iter().find(|&&(_, len)| len == n)?.0;

    // Identify principal functional group & bond unsaturation along best_chain
    let in_chain: HashMap<usize, usize> = best_chain.iter().enumerate().map(|(i, &a)| (a, i)).collect();
    let mut best_group: Option<(u8, usize, &str)> = None; // (priority, chain_idx, suffix)
    let mut dbl_locs: Vec<usize> = Vec::new();
    let mut trp_locs: Vec<usize> = Vec::new();

    for (idx, &ci) in best_chain.iter().enumerate() {
        if idx + 1 < n {
            let nxt = best_chain[idx + 1];
            if let Some(b) = mol.bonds.iter().find(|b| {
                (b.begin_atom == ci && b.end_atom == nxt) || (b.begin_atom == nxt && b.end_atom == ci)
            }) {
                if b.order == 2 { dbl_locs.push(idx); }
                if b.order == 3 { trp_locs.push(idx); }
            }
        }
        let ext_nbrs: Vec<(usize, i32)> = mol
            .neighbors(ci)
            .into_iter()
            .filter(|&(nb, _)| !in_chain.contains_key(&nb) && mol.atoms[nb].z != 1)
            .map(|(nb, bi)| (nb, mol.bonds[bi].order))
            .collect();
        let has_dbl_o = ext_nbrs.iter().any(|&(nb, ord)| mol.atoms[nb].z == 8 && ord == 2);
        let has_oh = ext_nbrs.iter().any(|&(nb, ord)| mol.atoms[nb].z == 8 && ord == 1 && mol.neighbors(nb).iter().filter(|&&(x, _)| mol.atoms[x].z != 1).count() == 1);
        let has_nh2 = ext_nbrs.iter().any(|&(nb, ord)| mol.atoms[nb].z == 7 && ord == 1);
        let has_cn = ext_nbrs.iter().any(|&(nb, ord)| mol.atoms[nb].z == 7 && ord == 3);

        let cand = if has_dbl_o && has_oh {
            Some((10, idx, "oic acid"))
        } else if has_dbl_o && has_nh2 {
            Some((9, idx, "amide"))
        } else if has_cn {
            Some((8, idx, "nitrile"))
        } else if has_dbl_o && (idx == 0 || idx == n - 1) {
            Some((7, idx, "al"))
        } else if has_dbl_o {
            Some((6, idx, "one"))
        } else if has_oh {
            Some((5, idx, "ol"))
        } else if has_nh2 {
            Some((4, idx, "amine"))
        } else {
            None
        };
        if let Some(c) = cand {
            if best_group.map(|b| c.0 > b.0).unwrap_or(true) {
                best_group = Some(c);
            }
        }
    }

    // Orient chain numbering to give lowest locant to principal functional group, unsaturation, or substituents!
    let reverse = if let Some((_, gidx, _)) = best_group {
        gidx > (n - 1) / 2
    } else if let Some(&d) = dbl_locs.first() {
        d >= n / 2
    } else if let Some(&t) = trp_locs.first() {
        t >= n / 2
    } else {
        let mut fwd_locs: Vec<usize> = Vec::new();
        let mut rev_locs: Vec<usize> = Vec::new();
        for (idx, &ci) in best_chain.iter().enumerate() {
            let sub_cnt = mol
                .neighbors(ci)
                .into_iter()
                .filter(|&(nb, _)| !in_chain.contains_key(&nb) && mol.atoms[nb].z != 1)
                .count();
            for _ in 0..sub_cnt {
                fwd_locs.push(idx + 1);
                rev_locs.push(n - idx);
            }
        }
        fwd_locs.sort_unstable();
        rev_locs.sort_unstable();
        rev_locs < fwd_locs
    };
    if reverse {
        best_chain.reverse();
    }
    let in_chain: HashMap<usize, usize> = best_chain.iter().enumerate().map(|(i, &a)| (a, i)).collect();
    let grp_pos = best_group.map(|(_, idx, suf)| (if reverse { n - 1 - idx } else { idx } + 1, suf));

    // Collect non-principal substituents
    let mut subs: Vec<(usize, String)> = Vec::new();
    for (idx, &ci) in best_chain.iter().enumerate() {
        let loc = idx + 1;
        for &(nb, bi) in &mol.neighbors(ci) {
            if in_chain.contains_key(&nb) || mol.atoms[nb].z == 1 {
                continue;
            }
            // Skip atoms belonging to the principal functional group on this carbon
            if let Some((gloc, suf)) = grp_pos {
                if gloc == loc {
                    let z = mol.atoms[nb].z;
                    let ord = mol.bonds[bi].order;
                    if (matches!(suf, "oic acid" | "al" | "one" | "ol") && z == 8)
                        || (matches!(suf, "amide") && (z == 8 || z == 7))
                        || (matches!(suf, "amine" | "nitrile") && z == 7 && (ord == 1 || ord == 3))
                    {
                        continue;
                    }
                }
            }
            subs.push((loc, classify_substituent(mol, nb, ci)));
        }
    }

    let infix = if !trp_locs.is_empty() {
        let l = if reverse { n - 1 - trp_locs[0] } else { trp_locs[0] + 1 };
        if n <= 3 { "yn".to_string() } else { format!("-{l}-yn") }
    } else if !dbl_locs.is_empty() {
        let l = if reverse { n - 1 - dbl_locs[0] } else { dbl_locs[0] + 1 };
        if n <= 3 { "en".to_string() } else { format!("-{l}-en") }
    } else {
        "an".to_string()
    };

    let parent = match grp_pos {
        None => format!("{root}{infix}e"),
        Some((loc, suf)) => match suf {
            "nitrile" => format!("{root}{infix}e{suf}"),
            "oic acid" | "al" | "amide" => format!("{root}{infix}{suf}"),
            _ => {
                if n <= 2 {
                    format!("{root}{infix}{suf}")
                } else {
                    format!("{root}{infix}-{loc}-{suf}")
                }
            }
        },
    };
    let pref = format_substituents(&subs);
    Some(format!("{pref}{parent}"))
}

fn format_substituents(subs: &[(usize, String)]) -> String {
    if subs.is_empty() {
        return String::new();
    }
    let mut grouped: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (loc, name) in subs {
        grouped.entry(name.clone()).or_default().push(*loc);
    }
    let mut parts = Vec::new();
    for (name, mut locs) in grouped {
        locs.sort_unstable();
        let mult = match locs.len() {
            1 => "",
            2 => "di",
            3 => "tri",
            4 => "tetra",
            _ => "poly",
        };
        let loc_s: Vec<String> = locs.iter().map(|l| l.to_string()).collect();
        parts.push(format!("{}-{mult}{name}", loc_s.join(",")));
    }
    parts.join("-")
}

fn dfs_longest_c_chain(
    mol: &ChemMolecule,
    u: usize,
    visited: &mut [bool],
    cur: &mut Vec<usize>,
    best: &mut Vec<usize>,
    budget: &mut usize,
) {
    if *budget == 0 {
        return;
    }
    *budget -= 1;
    if cur.len() > best.len() {
        *best = cur.clone();
    }
    if cur.len() >= 13 {
        return;
    }
    for &(v, _) in &mol.neighbors(u) {
        if !visited[v] && mol.atoms[v].z == 6 {
            visited[v] = true;
            cur.push(v);
            dfs_longest_c_chain(mol, v, visited, cur, best, budget);
            cur.pop();
            visited[v] = false;
        }
    }
}

/// Convenience helper for converting a molecule with simple bond/atom construction in tests.
pub fn make_simple_chain(n: usize) -> ChemMolecule {
    let mut m = ChemMolecule::default();
    for i in 0..n {
        let mut a = ChemAtom::new((i + 1) as u32, 6);
        a.x = (i as f64) * 14.4;
        m.atoms.push(a);
        if i > 0 {
            m.bonds.push(ChemBond {
                id: (n + i + 1) as u32,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iupac_name_to_smiles_and_structure_to_name() {
        let dict = NameDictionary::with_builtins();
        // Built-ins
        assert_eq!(dict.name_to_smiles("Aspirin").unwrap(), "CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(dict.name_to_smiles("Caffeine").unwrap(), "Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        // Systematic acyclic & cyclic IUPAC names
        assert_eq!(dict.name_to_smiles("2-chlorobutane").unwrap(), "CC(Cl)CC");
        assert_eq!(dict.name_to_smiles("3-pentanone").unwrap(), "CCC(=O)CC");
        assert_eq!(dict.name_to_smiles("4-nitrobenzoic acid").unwrap(), "c1(C(=O)O)ccc([N+](=O)[O-])cc1");
        assert_eq!(dict.name_to_smiles("cyclohexanol").unwrap(), "C1(O)CCCCC1");
        assert_eq!(dict.name_to_smiles("ethyl propanoate").unwrap(), "CCOC(=O)CC");
        // Structure -> Name
        let pentane = make_simple_chain(5);
        assert_eq!(molecule_to_name(&pentane, "CCCCC", &dict), "pentane");

        // Lowest-locant orientation on substituted alkane (chlorine on index 3 of 4-C chain => 2-chlorobutane, not 3-chlorobutane!)
        let mut chlorobutane = make_simple_chain(4);
        let mut cl = ChemAtom::new(99, 17);
        cl.x = 43.2;
        cl.y = 14.4;
        let cl_i = chlorobutane.add_atom(cl);
        chlorobutane.add_bond(99, 2, cl_i, 1, BondDisplay::Solid);
        assert_eq!(molecule_to_name(&chlorobutane, "CCC(Cl)C", &dict), "2-chlorobutane");

        // Substituted cycloalkane (methylcyclohexanol)
        let mut cyc = ChemMolecule::new();
        for i in 0..6 {
            let ang = (i as f64) * std::f64::consts::TAU / 6.0;
            let mut c = ChemAtom::new((i + 1) as u32, 6);
            c.x = 14.4 * ang.cos();
            c.y = 14.4 * ang.sin();
            cyc.add_atom(c);
        }
        for i in 0..6 {
            cyc.add_bond((10 + i) as u32, i, (i + 1) % 6, 1, BondDisplay::Solid);
        }
        let oh = cyc.add_atom(ChemAtom::new(20, 8));
        cyc.add_bond(21, 0, oh, 1, BondDisplay::Solid);
        let me = cyc.add_atom(ChemAtom::new(22, 6));
        cyc.add_bond(23, 5, me, 1, BondDisplay::Solid); // adjacent to OH on CCW side -> 2-methylcyclohexanol!
        assert_eq!(molecule_to_name(&cyc, "CC1CCCCC1O", &dict), "2-methylcyclohexanol");
        assert!(!is_fallback_chemical_name("2-methylcyclohexanol"));
        assert!(is_fallback_chemical_name("Compound (C17H19NO3)"));
    }

    #[test]
    fn test_bridged_cages_and_macrocycle_naming() {
        let dict = NameDictionary::with_builtins();
        // Adamantane (10 C, 12 single bonds, four deg-3 bridgeheads + six deg-2 bridges)
        let mut adam = ChemMolecule::new();
        for i in 0..10 {
            adam.add_atom(ChemAtom::new((i + 1) as u32, 6));
        }
        // Bridgeheads: 0, 1, 2, 3; CH2 bridges: 4=(0,1), 5=(0,2), 6=(0,3), 7=(1,2), 8=(1,3), 9=(2,3)
        let edges = [
            (0, 4), (4, 1), (0, 5), (5, 2), (0, 6), (6, 3),
            (1, 7), (7, 2), (1, 8), (8, 3), (2, 9), (9, 3),
        ];
        for (idx, &(u, v)) in edges.iter().enumerate() {
            adam.add_bond((20 + idx) as u32, u, v, 1, BondDisplay::Solid);
        }
        assert_eq!(molecule_to_name(&adam, "", &dict), "adamantane");

        // Norbornane (7 C, 8 single bonds, two deg-3 bridgeheads + five deg-2)
        let mut norb = ChemMolecule::new();
        for i in 0..7 {
            norb.add_atom(ChemAtom::new((i + 1) as u32, 6));
        }
        let nedges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (6, 3)];
        for (idx, &(u, v)) in nedges.iter().enumerate() {
            norb.add_bond((40 + idx) as u32, u, v, 1, BondDisplay::Solid);
        }
        assert_eq!(molecule_to_name(&norb, "", &dict), "norbornane (bicyclo[2.2.1]heptane)");

        // Cubane (8 C, 12 single bonds, all deg-3)
        let mut cub = ChemMolecule::new();
        for i in 0..8 {
            cub.add_atom(ChemAtom::new((i + 1) as u32, 6));
        }
        let cedges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        for (idx, &(u, v)) in cedges.iter().enumerate() {
            cub.add_bond((60 + idx) as u32, u, v, 1, BondDisplay::Solid);
        }
        assert_eq!(molecule_to_name(&cub, "", &dict), "cubane");
    }
}
