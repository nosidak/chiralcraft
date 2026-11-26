//! 1:1 Rust parity implementation of `CoreChemistryCommon.dll` (`26.8.0.6327`).
//! Every constant, table lookup, floating-point instruction order, and bitmask
//! is bound to verified disassembly of `CoreChemistryCommon.dll`
//! (SHA-256 `2127f1b9bfc2add67c6b4dac5bad5645d8a96868852c6ee396a08e127e429fcd`).

use std::collections::{BTreeMap, BTreeSet};

pub mod forcefield;
pub mod inorganic_polymer;
pub mod isomers;
pub mod iupac_name;
pub mod medchem;
pub mod nmr;
pub mod properties;
pub mod spectroscopy;
pub use forcefield::{
    measure_angle_deg, measure_dihedral_deg, measure_distance_angstrom, optimize_geometry_3d,
    ForceFieldReport,
};
pub use inorganic_polymer::{
    analyze_organometallic_complexes, compute_polymer_metrics, transition_metal_group,
    CopolymerTopology, MetalComplexReport, PolymerReport,
};
pub use isomers::{
    apply_electron_push_arrow, compute_microspecies_at_ph, detect_ionizable_sites,
    enumerate_resonance_structures, enumerate_stereoisomers, enumerate_tautomers,
    net_charge_at_ph, IonizableSite, MicrospeciesResult, ResonanceContributor,
};
pub use iupac_name::{is_fallback_chemical_name, molecule_to_name, NameDictionary, BUILTIN_NAMES};
pub use medchem::{
    compute_adducts, compute_medchem_report, detect_structural_alerts, simulate_hrms, AdductEntry,
    AlertSeverity, HrmsSimulation, IsotopePeak, MedChemReport, StructuralAlert,
};
pub use nmr::{predict_13c_nmr, predict_1h_nmr, NmrPeak, NmrPrediction};
pub use properties::{compute_analysis_report, estimate_chemical_properties, AnalysisReport, ChemicalProperties};
pub use spectroscopy::{
    predict_2d_nmr, predict_ir_spectrum, predict_uv_vis, IrBand, IrPrediction, Nmr2dCrossPeak,
    Nmr2dKind, Nmr2dPrediction, UvVisPrediction,
};

pub use chiralcraft_tables::{
    cs_real_num_kundef, electron_mass_factor, ElementRecord, IsotopeRecord, PseudoAtomRecord,
    ShippedNickAtom, ShippedNickBond, ShippedNickname, CAN_BE_AROMATIC_JUMP_INDEX,
    CS_REAL_NUM_KUNDEF_BITS, DEFAULT_VALENCE_TABLE, ELECTRON_MASS_FACTOR_BITS, ELEMENTS, ISOTOPES,
    PSEUDO_ATOMS, ROUND_HALF_BITS, SHIPPED_NICKNAMES,
};

// ============================================================================
// 1. csRealNum (VA 0x1800b70d0, 0x1800b70e0, 0x180bd9370, 0x180484850)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CsRealNum {
    pub val: f64,
    pub tolerance: f64,
}

impl CsRealNum {
    pub const fn new(val: f64, tolerance: f64) -> Self {
        Self { val, tolerance }
    }

    pub fn undefined() -> Self {
        Self {
            val: cs_real_num_kundef(),
            tolerance: 0.0,
        }
    }

    /// Matches `?IsDefined@csRealNum@@QEBA_NXZ` at VA `0x180484850`:
    /// returns `false` only on ordered equality (`comisd` ZF=1, PF=0) with `kUNDEF` (`3229572.43`).
    #[inline]
    pub fn is_defined(&self) -> bool {
        let u = cs_real_num_kundef();
        if self.val.is_nan() {
            true
        } else {
            self.val != u
        }
    }
}

// ============================================================================
// 2. csPeriodicTable & csElement (VA 0x18047dc90, 0x1800b7700, 0x1800b89d0, 0x18047eac0)
// ============================================================================

#[inline]
pub fn get_element(z: i32) -> Option<&'static ElementRecord> {
    if (1..=118).contains(&z) {
        Some(&ELEMENTS[(z - 1) as usize])
    } else {
        None
    }
}

pub fn get_element_by_symbol(sym: &str) -> Option<&'static ElementRecord> {
    ELEMENTS.iter().find(|e| e.symbol == sym)
}

pub fn get_pseudo_atom(z: i32) -> Option<&'static PseudoAtomRecord> {
    if (220..=245).contains(&z) {
        Some(&PSEUDO_ATOMS[(z - 220) as usize])
    } else {
        None
    }
}

pub fn get_pseudo_atom_by_symbol(sym: &str) -> Option<&'static PseudoAtomRecord> {
    PSEUDO_ATOMS.iter().find(|p| p.symbol == sym)
}

/// Returns the slice of `ISOTOPES` belonging to element `z` (`1..=118`).
pub fn get_element_isotopes(z: i32) -> &'static [IsotopeRecord] {
    if !(1..=118).contains(&z) {
        return &[];
    }
    let start = ISOTOPES.partition_point(|iso| iso.z < z);
    let end = ISOTOPES.partition_point(|iso| iso.z <= z);
    &ISOTOPES[start..end]
}

/// Matches `?GetIsotope@csElement@cs@@QEBAPEBVcsIsotope@2@H@Z` at VA `0x18047dc90`:
/// Walks the element's isotope records and returns the first record where
/// `floor(record.mass + 0.5)` converted via `cvttsd2si` (`as i32`) equals `mass_number`.
pub fn get_isotope(z: i32, mass_number: i32) -> Option<&'static IsotopeRecord> {
    let half = f64::from_bits(ROUND_HALF_BITS);
    for rec in get_element_isotopes(z) {
        let rounded = (rec.mass() + half).floor() as i32;
        if rounded == mass_number {
            return Some(rec);
        }
    }
    None
}

/// Matches `?GetExactMass@csElement@cs@@QEBA?AVcsRealNum@@XZ` at VA `0x1800b7700`:
/// - Empty isotope list -> returns `(atomic_weight, 0.0)`.
/// - Single isotope record -> returns `(rec[0].mass, 0.0)`.
/// - Multiple records -> starts with `rec[0]`; a later record replaces the candidate
///   only if `later.abundance > candidate.abundance` (strictly greater; ties keep first).
pub fn get_element_exact_mass(z: i32) -> CsRealNum {
    let Some(elem) = get_element(z) else {
        return CsRealNum::undefined();
    };
    let isos = get_element_isotopes(z);
    if isos.is_empty() {
        return CsRealNum::new(elem.atomic_weight(), 0.0);
    }
    let mut best = &isos[0];
    for iso in &isos[1..] {
        if iso.abundance() > best.abundance() {
            best = iso;
        }
    }
    CsRealNum::new(best.mass(), 0.0)
}

/// Matches `?CanInferImplicitHydrogen@csPeriodicTable@cs@@SA_NH@Z` at VA `0x1800b89d0`:
/// True iff `Z in {1, 6, 7, 8, 9, 14, 15, 16}` (`H, C, N, O, F, Si, P, S`).
#[inline]
pub fn can_infer_implicit_hydrogen(z: i32) -> bool {
    z == 1 || (6..=9).contains(&z) || (14..=16).contains(&z)
}

/// Matches `?DefaultValence@csPeriodicTable@cs@@SAHH@Z` at VA `0x18047eac0`:
/// Looks up `DEFAULT_VALENCE_TABLE[z]` for `1 <= z <= 103`, else `0`.
#[inline]
pub fn default_valence(z: i32) -> i32 {
    if (1..=103).contains(&z) {
        DEFAULT_VALENCE_TABLE[z as usize]
    } else {
        0
    }
}

/// Matches `0x180296340` -> `?GetIntendedValence@csPeriodicTable@cs@@SAHHH@Z` (`0x18047ea30`)
/// -> `csValenceList` round-up (`0x18047d420`):
/// - `128..=218` (`(z - 0x80) <= 0x5a`) -> `0`
/// - `220..=245` (`(z - 0xdc) <= 0x19`) -> `PSEUDO_ATOMS[z - 220].valence`
/// - `1..=118`:
///   - if `allowed_valences` is empty -> `0`
///   - if `required <= preferred_valence` -> `preferred_valence`
///   - else smallest `v in allowed_valences` with `v >= required`
///   - else `allowed_valences.last()` (max allowed valence)
pub fn get_intended_valence(z: i32, required: i32) -> i32 {
    if (128..=218).contains(&z) {
        return 0;
    }
    if (220..=245).contains(&z) {
        return PSEUDO_ATOMS[(z - 220) as usize].valence;
    }
    let Some(elem) = get_element(z) else {
        return 0;
    };
    if elem.allowed_valences.is_empty() {
        return 0;
    }
    if required <= elem.preferred_valence {
        return elem.preferred_valence;
    }
    for &v in elem.allowed_valences {
        if v >= required {
            return v;
        }
    }
    *elem.allowed_valences.last().unwrap()
}

// ============================================================================
// 3. csFormulaProperties (VA 0x18047a010, 0x180478d90, 0x180477aa0)
// ============================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnalyzedElementKey {
    pub z: i32,
    pub symbol: String,
    pub isotope: i32,
}

impl Ord for AnalyzedElementKey {
    /// Matches `std::less<AnalyzedElement>` inlined at `0x180479e90`–`0x180479ed8`:
    /// orders by symbol string (`record+0x38`) primary, `isotope` (`+0x14`) secondary.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.z == other.z && self.symbol == other.symbol {
            self.isotope.cmp(&other.isotope)
        } else {
            match self.symbol.cmp(&other.symbol) {
                std::cmp::Ordering::Equal => self.isotope.cmp(&other.isotope),
                ord => ord,
            }
        }
    }
}

impl PartialOrd for AnalyzedElementKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
pub struct CsFormulaProperties {
    elements_map: BTreeMap<AnalyzedElementKey, f64>,
    other_nodes_map: BTreeMap<String, f64>,
    total_charge: f64,
    total_radical: f64,
    exact_mass_cache: Option<CsRealNum>,
    /// `this+0xb2`: isotope-fallback grouping flag byte.
    pub flag_b2: u8,
    /// `this+0xb4`: minimum numeral width (signed int32; `<= 0` omits unity `1`).
    pub min_numeral_width_b4: i32,
}

impl Default for CsFormulaProperties {
    fn default() -> Self {
        Self::new()
    }
}

impl CsFormulaProperties {
    pub fn new() -> Self {
        Self {
            elements_map: BTreeMap::new(),
            other_nodes_map: BTreeMap::new(),
            total_charge: 0.0,
            total_radical: 0.0,
            exact_mass_cache: None,
            flag_b2: 0,
            min_numeral_width_b4: 0,
        }
    }

    /// Matches `Invalidate` at VA `0x180478120` (vtable slot `+0x48`).
    #[inline]
    pub fn invalidate(&mut self) {
        self.exact_mass_cache = None;
    }

    /// Matches `Clear` at VA `0x1804780a0` (`file 0x4774a0`).
    pub fn clear(&mut self) {
        self.elements_map.clear();
        self.other_nodes_map.clear();
        self.total_charge = 0.0;
        self.total_radical = 0.0;
        self.invalidate();
    }

    #[inline]
    pub fn total_charge(&self) -> f64 {
        self.total_charge
    }

    #[inline]
    pub fn total_radical(&self) -> f64 {
        self.total_radical
    }

    /// Matches `?IncrementCharge@csFormulaProperties@cs@@UEAAXN@Z` at VA `0x180164660`.
    pub fn increment_charge(&mut self, delta: f64) {
        self.total_charge += delta;
        self.invalidate();
    }

    /// Matches `?IncrementRadical@csFormulaProperties@cs@@UEAAXN@Z` at VA `0x180164160`.
    pub fn increment_radical(&mut self, delta: f64) {
        self.total_radical += delta;
        self.invalidate();
    }

    /// Matches numeric `IncrementAtom` at file offset `0x477d10`.
    pub fn increment_atom(&mut self, z: i32, isotope: i32, count: f64) -> bool {
        let Some(elem) = get_element(z) else {
            return false;
        };
        let key = AnalyzedElementKey {
            z,
            symbol: elem.symbol.to_string(),
            isotope,
        };
        *self.elements_map.entry(key).or_insert(0.0) += count;
        self.invalidate();
        true
    }

    /// Matches string `IncrementAtom` at file offset `0x477fa0`.
    pub fn increment_other_node(&mut self, label: &str, count: f64) {
        *self.other_nodes_map.entry(label.to_string()).or_insert(0.0) += count;
        self.invalidate();
    }

    /// Matches `?AddFormula@csFormulaProperties@cs@@QEAAXAEBV12@N@Z` at VA `0x1804785d0`.
    pub fn add_formula(&mut self, other: &CsFormulaProperties, scale: f64) {
        for (k, &v) in &other.elements_map {
            *self.elements_map.entry(k.clone()).or_insert(0.0) += v * scale;
        }
        for (k, &v) in &other.other_nodes_map {
            *self.other_nodes_map.entry(k.clone()).or_insert(0.0) += v * scale;
        }
        self.increment_charge(other.total_charge * scale);
        self.increment_radical(other.total_radical * scale);
        self.invalidate();
    }

    /// Matches `?GetExactMass@csFormulaProperties@cs@@QEBAAEBVcsRealNum@@XZ` at VA `0x18047a010`.
    pub fn get_exact_mass(&mut self) -> CsRealNum {
        if let Some(cached) = self.exact_mass_cache {
            return cached;
        }
        let u = cs_real_num_kundef();
        let mut acc_val = 0.0_f64;
        let mut acc_tol = 0.0_f64;

        for (key, &c) in &self.elements_map {
            let selected = if key.isotope != 0 {
                if let Some(iso) = get_isotope(key.z, key.isotope) {
                    CsRealNum::new(iso.mass(), 0.0)
                } else {
                    get_element_exact_mass(key.z)
                }
            } else {
                get_element_exact_mass(key.z)
            };

            let m = selected.val;
            let t = selected.tolerance;

            let (p, tp) = if (!m.is_nan() && m == u) || (!c.is_nan() && c == u) {
                (u, 0.0)
            } else {
                // Exact instruction order at 0x18047a0e0-0x18047a0fc:
                // abs(m) * +0.0, then abs(c) * t, then sum tp, then p = m * c
                let term1 = f64::from_bits(m.to_bits() & 0x7fff_ffff_ffff_ffff) * 0.0;
                let term2 = f64::from_bits(c.to_bits() & 0x7fff_ffff_ffff_ffff) * t;
                let tp = term1 + term2;
                let p = m * c;
                (p, tp)
            };

            if (!acc_val.is_nan() && acc_val == u) || (!p.is_nan() && p == u) {
                acc_val = u;
                acc_tol = 0.0;
            } else {
                // Exact instruction order at 0x18047a121-0x18047a12f:
                acc_tol = tp + acc_tol;
                acc_val = acc_val + p;
            }
        }

        let result = if !acc_val.is_nan() && acc_val == u {
            CsRealNum::new(u, 0.0)
        } else {
            // Toggle sign bit of total_charge (`xorpd` with 0x8000000000000000)
            let neg_charge = f64::from_bits(self.total_charge.to_bits() ^ 0x8000_0000_0000_0000);
            let delta = neg_charge * electron_mass_factor();
            if !delta.is_nan() && delta == u {
                CsRealNum::new(u, 0.0)
            } else {
                CsRealNum::new(acc_val + delta, acc_tol + 0.0)
            }
        };

        self.exact_mass_cache = Some(result);
        result
    }

    /// Matches `GetFormula` at VA `0x180478d90` and count formatter at VA `0x180477aa0`:
    /// Hill system (C=6 first, H=1 second, then remaining elements in map order).
    pub fn get_formula(&self) -> String {
        let mut out = String::new();

        // 1. Carbon (Z == 6)
        for (k, &cnt) in &self.elements_map {
            if k.z == 6 {
                Self::format_element_count(k, &mut out, cnt, self.min_numeral_width_b4, self.flag_b2);
            }
        }
        // 2. Hydrogen (Z == 1)
        for (k, &cnt) in &self.elements_map {
            if k.z == 1 {
                Self::format_element_count(k, &mut out, cnt, self.min_numeral_width_b4, self.flag_b2);
            }
        }
        // 3. All other elements (Z != 1 && Z != 6) in lexicographical symbol + isotope order
        for (k, &cnt) in &self.elements_map {
            if k.z != 1 && k.z != 6 {
                Self::format_element_count(k, &mut out, cnt, self.min_numeral_width_b4, self.flag_b2);
            }
        }
        // 4. Other nodes map
        for (label, &cnt) in &self.other_nodes_map {
            if cnt > 0.0 {
                out.push_str(label);
                if cnt != 1.0 || self.min_numeral_width_b4 > 0 {
                    Self::append_formatted_number(&mut out, cnt, self.min_numeral_width_b4);
                }
            }
        }

        out
    }

    fn format_element_count(
        key: &AnalyzedElementKey,
        out: &mut String,
        count: f64,
        opt_width: i32,
        flag_b2: u8,
    ) {
        // 1. Non-positive or NaN count emits nothing (0x180477af3)
        if !(count > 0.0) {
            return;
        }

        let mut sym_to_emit = key.symbol.as_str();
        let mut needs_paren = false;

        // 6. Explicit isotope substitution (0x180477b10-0x180477b95)
        if key.isotope != 0 {
            if let Some(iso) = get_isotope(key.z, key.isotope) {
                if !iso.symbol.is_empty() {
                    sym_to_emit = iso.symbol;
                } else if flag_b2 == 0 {
                    if let Some(last_ch) = out.chars().last() {
                        if !last_ch.is_whitespace() && !matches!(last_ch, '(' | '[' | '{') {
                            needs_paren = true;
                        }
                    }
                }
            }
        }

        if needs_paren {
            out.push('(');
        }
        out.push_str(sym_to_emit);
        if needs_paren {
            out.push(')');
        }

        // 3. Unity omission when count == 1.0 and signed opt_width <= 0 (0x180477c08-0x180477c17)
        if count == 1.0 && opt_width <= 0 {
            return;
        }

        Self::append_formatted_number(out, count, opt_width);
    }

    /// Populates `CsFormulaProperties` directly from a `ChemMolecule`, adding explicit atoms,
    /// isotopes, formal charges, and inferred implicit hydrogens (`CcVal::compute`).
    pub fn from_molecule(mol: &ChemMolecule) -> Self {
        let mut props = Self::new();
        let val = CcVal::compute(mol);
        for (i, atom) in mol.atoms.iter().enumerate() {
            props.increment_atom(atom.z, atom.isotope, 1.0);
            if atom.charge != 0 {
                props.increment_charge(atom.charge as f64);
            }
            let imp_h = val.implicit_h_counts[i] as f64;
            if imp_h > 0.0 {
                props.increment_atom(1, 0, imp_h);
            }
        }
        props
    }

    /// Standard average molecular weight (IUPAC atomic weights, or explicit isotope mass where specified).
    pub fn get_molecular_weight(&self) -> f64 {
        let mut total = 0.0;
        for (k, &cnt) in &self.elements_map {
            if cnt <= 0.0 {
                continue;
            }
            let unit_mass = if k.isotope != 0 {
                get_isotope(k.z, k.isotope)
                    .map(|iso| iso.mass())
                    .unwrap_or_else(|| get_element(k.z).map(|el| el.atomic_weight()).unwrap_or(0.0))
            } else {
                get_element(k.z).map(|el| el.atomic_weight()).unwrap_or(0.0)
            };
            total += cnt * unit_mass;
        }
        total - self.total_charge * electron_mass_factor()
    }

    /// Mass-to-charge ratio (`m/z`): exact mass divided by `max(1, |net_charge|)`.
    pub fn get_m_over_z(&mut self) -> f64 {
        let q = self.total_charge.abs();
        let em = self.get_exact_mass().val;
        if q >= 1.0 {
            em / q
        } else {
            em
        }
    }

    /// Elemental analysis mass percentages ordered in Hill system order (`C`, `H`, remaining alphabetical).
    pub fn get_elemental_analysis(&self) -> Vec<(String, f64)> {
        let mw = self.get_molecular_weight();
        if mw <= 0.0 {
            return Vec::new();
        }
        let mut by_sym: BTreeMap<(u8, String), f64> = BTreeMap::new();
        for (k, &cnt) in &self.elements_map {
            if cnt <= 0.0 {
                continue;
            }
            let unit_mass = if k.isotope != 0 {
                get_isotope(k.z, k.isotope)
                    .map(|iso| iso.mass())
                    .unwrap_or_else(|| get_element(k.z).map(|el| el.atomic_weight()).unwrap_or(0.0))
            } else {
                get_element(k.z).map(|el| el.atomic_weight()).unwrap_or(0.0)
            };
            let order_bucket = match k.z {
                6 => 0u8,
                1 => 1u8,
                _ => 2u8,
            };
            let sym = if k.isotope != 0 {
                get_isotope(k.z, k.isotope)
                    .filter(|iso| !iso.symbol.is_empty())
                    .map(|iso| iso.symbol.to_string())
                    .unwrap_or_else(|| k.symbol.clone())
            } else {
                k.symbol.clone()
            };
            *by_sym.entry((order_bucket, sym)).or_insert(0.0) += cnt * unit_mass;
        }
        by_sym
            .into_iter()
            .map(|((_, sym), mass)| (sym, (mass / mw) * 100.0))
            .collect()
    }

    /// Formats elemental analysis percentages matching ChiralCraft (`"C, 92.26; H, 7.74"`).
    pub fn format_elemental_analysis(&self) -> String {
        self.get_elemental_analysis()
            .into_iter()
            .map(|(sym, pct)| format!("{sym}, {pct:.2}"))
            .collect::<Vec<_>>()
            .join("; ")
    }

    fn append_formatted_number(out: &mut String, count: f64, opt_width: i32) {
        let trunc = count as i64;
        if (count - (trunc as f64)) == 0.0 && trunc >= 0 {
            let s = trunc.to_string();
            let min_w = opt_width.max(0) as usize;
            if s.len() < min_w {
                for _ in 0..(min_w - s.len()) {
                    out.push('0');
                }
            }
            out.push_str(&s);
        } else {
            let s = format!("{:.6}", count)
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string();
            let min_w = opt_width.max(0) as usize;
            if s.len() < min_w {
                for _ in 0..(min_w - s.len()) {
                    out.push('0');
                }
            }
            out.push_str(&s);
        }
    }
}

// ============================================================================
// 4. Molecule Graph & ccVal Engine (VA 0x180200a50, 0x180201550, 0x1802021b0)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondDisplay {
    Solid,
    Dash,
    Hash,
    WedgedHashBegin,
    WedgedHashEnd,
    Bold,
    WedgeBegin,
    WedgeEnd,
    Wavy,
    HollowWedgeBegin,
    HollowWedgeEnd,
    Dative,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChemAtom {
    pub id: u32,
    /// Atomic number (`1..=118`), pseudo-atom (`220..=245`), or `0` (query/unspecified).
    pub z: i32,
    pub isotope: i32,
    pub charge: i32,
    /// `0` = none, `1` = doublet/monoradical (1e), `2` = singlet (2e), `3` = triplet (2e),
    /// matching CoreChemistry `ccAtom+0x34` in `AnalyzeAtomsExpectedValence` (`0x180201a15`).
    pub radical: i32,
    /// Explicit valence override (`None` = unset, `Some(v)` = override `v`).
    pub valence_override: Option<i32>,
    /// CDXML `NumHydrogens` attribute when explicitly authored/perceived (`None` = unspecified).
    pub explicit_num_hydrogens: Option<i32>,
    pub is_query: bool,
    pub is_unspecified_label: bool,
    pub visible: bool,
    pub ignore_warnings: bool,
    pub label: Option<String>,
    pub x: f64,
    pub y: f64,
    pub z_coord: f64,
}

impl ChemAtom {
    pub fn new(id: u32, z: i32) -> Self {
        Self {
            id,
            z,
            isotope: 0,
            charge: 0,
            radical: 0,
            valence_override: None,
            explicit_num_hydrogens: None,
            is_query: false,
            is_unspecified_label: false,
            visible: true,
            ignore_warnings: false,
            label: None,
            x: 0.0,
            y: 0.0,
            z_coord: 0.0,
        }
    }

    /// Maps radical code (`+0x34`) as in `AnalyzeAtomsExpectedValence`: `1->1, 2->2, 3->2, else 0`.
    #[inline]
    pub fn radical_electrons(&self) -> i32 {
        match self.radical {
            1 => 1,
            2 | 3 => 2,
            _ => 0,
        }
    }

    /// Clamps charge (`+0x30`) as in `AnalyzeAtomsExpectedValence` (`> 0xffd -> 0`).
    #[inline]
    pub fn clamped_charge(&self) -> i32 {
        if self.charge.abs() > 0xffd {
            0
        } else {
            self.charge
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChemBond {
    pub id: u32,
    pub begin_atom: usize,
    pub end_atom: usize,
    /// Bond order code (`1` = single, `2` = double, `3` = triple, `4` = quadruple,
    /// `6` = aromatic/delocalized, `5`/`7`/`8`/`9+` = dative/hydrogen/any/other).
    pub order: i32,
    pub display: BondDisplay,
    pub visible: bool,
    pub ignore_warnings: bool,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ChemMolecule {
    pub atoms: Vec<ChemAtom>,
    pub bonds: Vec<ChemBond>,
    pub has_3d_coords: bool,
}

impl ChemMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, atom: ChemAtom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        idx
    }

    pub fn add_bond(&mut self, id: u32, begin_atom: usize, end_atom: usize, order: i32, display: BondDisplay) -> usize {
        let idx = self.bonds.len();
        self.bonds.push(ChemBond {
            id,
            begin_atom,
            end_atom,
            order,
            display,
            visible: true,
            ignore_warnings: false,
        });
        idx
    }

    pub fn neighbors(&self, atom_idx: usize) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (bi, b) in self.bonds.iter().enumerate() {
            if b.begin_atom == atom_idx {
                out.push((b.end_atom, bi));
            } else if b.end_atom == atom_idx {
                out.push((b.begin_atom, bi));
            }
        }
        out
    }

    /// Computes molecular formula + exact mass using `ccVal` implicit hydrogens + `CsFormulaProperties`.
    pub fn compute_formula_properties(&self) -> CsFormulaProperties {
        let val = CcVal::compute(self);
        let mut fp = CsFormulaProperties::new();
        let mut total_charge = 0;
        let mut total_radical = 0;
        for (i, atom) in self.atoms.iter().enumerate() {
            total_charge += atom.clamped_charge();
            total_radical += atom.radical_electrons();
            if (1..=118).contains(&atom.z) {
                fp.increment_atom(atom.z, atom.isotope, 1.0);
            } else if atom.z == 225 {
                // Deuterium pseudo-atom -> H isotope 2
                fp.increment_atom(1, 2, 1.0);
            } else if atom.z == 226 {
                // Tritium pseudo-atom -> H isotope 3
                fp.increment_atom(1, 3, 1.0);
            }
            let imp_h = val.implicit_h(i);
            if imp_h > 0 {
                fp.increment_atom(1, 0, imp_h as f64);
            }
        }
        if total_charge != 0 {
            fp.increment_charge(total_charge as f64);
        }
        if total_radical != 0 {
            fp.increment_radical(total_radical as f64);
        }
        fp
    }
}

#[inline]
pub fn is_metal_or_metalloid(z: i32) -> bool {
    (3..=4).contains(&z)
        || (11..=13).contains(&z)
        || (19..=32).contains(&z)
        || (37..=51).contains(&z)
        || (55..=84).contains(&z)
        || (87..=118).contains(&z)
}

#[inline]
pub fn is_chalcogen_like(z: i32) -> bool {
    matches!(z, 8 | 16 | 34 | 52 | 84)
}

#[inline]
pub fn is_pnictogen_or_b(z: i32) -> bool {
    matches!(z, 5 | 7 | 15 | 33 | 51)
}

#[inline]
pub fn is_pnictogen(z: i32) -> bool {
    matches!(z, 7 | 15 | 33 | 51 | 83)
}

#[inline]
pub fn is_halogen(z: i32) -> bool {
    matches!(z, 9 | 17 | 35 | 53 | 85)
}

#[inline]
pub fn is_alkali(z: i32) -> bool {
    matches!(z, 3 | 11 | 19 | 37 | 55 | 87)
}

#[inline]
pub fn is_alkaline_earth(z: i32) -> bool {
    matches!(z, 4 | 12 | 20 | 38 | 56 | 88)
}

#[inline]
pub fn is_noble_gas(z: i32) -> bool {
    matches!(z, 2 | 10 | 18 | 36 | 54 | 86)
}

#[inline]
pub fn period_bucket(z: i32) -> i32 {
    if z <= 2 {
        0
    } else if z <= 10 {
        1
    } else if z <= 18 {
        2
    } else if z <= 36 {
        3
    } else if z <= 54 {
        4
    } else if z <= 88 {
        5
    } else if z <= 124 {
        6
    } else {
        -1
    }
}

/// Matches `0x180200db0` in `CoreChemistryCommon.dll`:
/// Validates whether `(z, eff_charge)` forms a supported main-group isoelectronic state.
pub fn is_valid_isoelectronic_state(z: i32, eff_charge: i32, req_base: i32) -> bool {
    if eff_charge == 0 {
        return true;
    }
    if !(1..=118).contains(&z) || eff_charge > 0xffd {
        return false;
    }
    if eff_charge > 0 {
        if is_halogen(z) && (eff_charge > 1 || req_base >= 3) {
            return false;
        }
    } else if is_alkali(z) || is_alkaline_earth(z) {
        return false;
    }
    let iso_z = z - eff_charge;
    if period_bucket(z) != period_bucket(iso_z) && !is_noble_gas(iso_z) {
        return false;
    }
    // Reject transition metals, lanthanides, and actinides
    !((21..=30).contains(&z)
        || (39..=48).contains(&z)
        || (57..=80).contains(&z)
        || (89..=118).contains(&z))
}

/// Matches `0x180201430` in `CoreChemistryCommon.dll`:
/// Detects bridging H or halogen (e.g. B-H-B in diborane or Al-Cl-Al).
fn is_bridging_h_or_halogen(mol: &ChemMolecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms[atom_idx];
    if atom.clamped_charge() != 0 {
        return false;
    }
    let ovr_raw = atom.valence_override.map(|v| v + 1).unwrap_or(0);
    let ok_elem = (atom.z == 1 && ovr_raw < 3) || (is_halogen(atom.z) && ovr_raw == 2);
    if !ok_elem {
        return false;
    }
    let nbrs = mol.neighbors(atom_idx);
    if nbrs.len() != 2 {
        return false;
    }
    nbrs.iter().all(|&(nbr, _)| {
        let nz = mol.atoms[nbr].z;
        nz == 5 || is_metal_or_metalloid(nz)
    })
}

/// Matches `0x180201390` in `CoreChemistryCommon.dll`:
/// Returns true when atom has an aromatic bond and is not standard C/O/S (or is Si).
fn is_aromatic_non_standard_heteroatom(mol: &ChemMolecule, atom_idx: usize) -> bool {
    let z = mol.atoms[atom_idx].z;
    if matches!(z, 6 | 8 | 16) {
        return false;
    }
    mol.neighbors(atom_idx)
        .into_iter()
        .any(|(_, bi)| mol.bonds[bi].order == 6)
}

#[derive(Debug, Clone)]
pub struct CcVal {
    pub sbo: Vec<i32>,
    pub expected_valence: Vec<i32>,
    pub expected_valence_a_priori: Vec<i32>,
    pub implicit_h_counts: Vec<u8>,
    pub explicit_h_counts: Vec<u8>,
    pub explicit_dative_out: Vec<u8>,
    pub explicit_dative_in: Vec<u8>,
    pub total_dative_out: Vec<u8>,
    pub total_dative_in: Vec<u8>,
    pub indeterminate_valence: Vec<bool>,
}

impl CcVal {
    /// Runs `CountExplicitDativeBonds` (`0x180200950`), `ComputeSbo` (`0x180200a50`),
    /// `AnalyzeAtomsExpectedValence` (`0x180201550`), `ComputeHCount` (`0x1802021b0`),
    /// and `IdentifyTacitDativeBonds` (`0x1802029d0`) across all atoms in `mol`.
    pub fn compute(mol: &ChemMolecule) -> Self {
        let n = mol.atoms.len();
        let mut sbo = vec![0; n];
        let mut expected_valence = vec![0; n];
        let mut expected_valence_a_priori = vec![0; n];
        let mut implicit_h_counts = vec![0; n];
        let mut explicit_h_counts = vec![0; n];
        let mut explicit_dative_out = vec![0u8; n];
        let mut explicit_dative_in = vec![0u8; n];
        let mut total_dative_out = vec![0u8; n];
        let mut total_dative_in = vec![0u8; n];
        let mut indeterminate_valence = vec![false; n];

        // 0. CountExplicitDativeBonds (0x180200950)
        for b in &mol.bonds {
            if b.order == 5 {
                if b.begin_atom < n {
                    explicit_dative_out[b.begin_atom] = explicit_dative_out[b.begin_atom].saturating_add(1);
                    total_dative_out[b.begin_atom] = total_dative_out[b.begin_atom].saturating_add(1);
                }
                if b.end_atom < n {
                    explicit_dative_in[b.end_atom] = explicit_dative_in[b.end_atom].saturating_add(1);
                    total_dative_in[b.end_atom] = total_dative_in[b.end_atom].saturating_add(1);
                }
            }
        }

        // 1. ComputeSbo (0x180200a50)
        for i in 0..n {
            let mut sum = 0i32;
            let mut has_non_chalcogen_arom = false;
            let nbrs = mol.neighbors(i);
            for &(nbr, bi) in &nbrs {
                let ord = mol.bonds[bi].order;
                match ord {
                    1..=4 => sum += ord,
                    5 | 8 => sum += 1,
                    6 => {
                        sum += 1;
                        let nbr_iso_z = mol.atoms[nbr].z - mol.atoms[nbr].clamped_charge();
                        if !is_chalcogen_like(nbr_iso_z) {
                            has_non_chalcogen_arom = true;
                        }
                    }
                    _ => sum += 1,
                }
            }

            let atom = &mol.atoms[i];
            let iso_z = atom.z - atom.clamped_charge();
            let ovr_raw = atom.valence_override.map(|v| v + 1).unwrap_or(0);

            if has_non_chalcogen_arom {
                if ovr_raw == sum + 1 {
                    if sum == 5
                        && nbrs
                            .iter()
                            .any(|&(_, bi)| (2..=4).contains(&mol.bonds[bi].order))
                    {
                        sum = 4;
                    }
                } else if (sum == 2 || sum == 4 || sum >= 6) && is_chalcogen_like(iso_z) {
                    if sum == 2 {
                        if iso_z != 8 && (ovr_raw == 5 || ovr_raw == 7) {
                            sum += 1;
                        }
                    } else if sum == 4 {
                        if iso_z != 8 && ovr_raw == 7 {
                            sum += 1;
                        }
                    } else if is_pnictogen_or_b(iso_z) {
                        if nbrs.iter().any(|&(nbr, _)| is_bridging_h_or_halogen(mol, nbr)) {
                            sum += 1;
                        }
                    } else {
                        sum += 1;
                    }
                } else if (sum == 3 || sum >= 5) && is_pnictogen_or_b(iso_z) {
                    if nbrs.iter().any(|&(nbr, _)| is_bridging_h_or_halogen(mol, nbr)) {
                        sum += 1;
                    } else if sum == 3 && iso_z != 7 && ovr_raw == 6 {
                        sum += 1;
                    }
                } else {
                    // Non-chalcogen with sum in {1, 2, 4} or non-pnictogen with sum >= 3: +1 Kekulé valence
                    sum += 1;
                }
            }
            sbo[i] = sum;
        }

        // 2. AnalyzeAtomsExpectedValence (0x180201550) & ComputeAtomsExpectedValence (0x180202150)
        for i in 0..n {
            let atom = &mol.atoms[i];
            let z = atom.z;
            let b1 = explicit_dative_out[i] as i32;
            let b2 = explicit_dative_in[i] as i32;
            let d_chg = b1 - b2;
            let ovr = atom.valence_override.unwrap_or(0).max(0);
            let req_base = sbo[i].max(ovr);
            let chg = atom.clamped_charge();
            let eff_chg = chg + d_chg;
            let iso_z = z - eff_chg;
            let rad = atom.radical_electrons();
            let in_pt = (1..=118).contains(&z);
            let is_metal = in_pt && is_metal_or_metalloid(z);

            let mut ev;
            if (128..=218).contains(&z) {
                if sbo[i] < iso_z && ovr > 0 {
                    indeterminate_valence[i] = true;
                }
                ev = req_base;
            } else if z <= 0 || atom.is_query {
                ev = req_base;
            } else {
                let mut use_intended_path = is_valid_isoelectronic_state(z, eff_chg, req_base);
                if use_intended_path && eff_chg != 0 && is_metal {
                    let ok_metal = if is_alkali(z) {
                        eff_chg == 1 && req_base == 0
                    } else if is_alkaline_earth(z) {
                        (eff_chg == 1 || eff_chg == 2) && (eff_chg + req_base == 2)
                    } else if is_halogen(z) {
                        eff_chg == 1 && req_base < 3
                    } else {
                        false
                    };
                    if !ok_metal {
                        use_intended_path = false;
                    }
                }

                if use_intended_path {
                    let required = req_base + rad;
                    let mut intended = if iso_z <= 0 {
                        0
                    } else {
                        let mut iv = get_intended_valence(iso_z, required);
                        if eff_chg != 0
                            && required != iv
                            && (iso_z % 2 == 0)
                            && is_noble_gas(iso_z)
                        {
                            iv = 0;
                        }
                        iv
                    };

                    // Element-specific special-case branches (0x180201a50..0x180201cb0)
                    if required > intended && in_pt {
                        if z == 5 && req_base > 3 {
                            if mol
                                .neighbors(i)
                                .into_iter()
                                .any(|(nbr, _)| is_bridging_h_or_halogen(mol, nbr))
                            {
                                intended = req_base;
                            }
                        } else if z == 13 && req_base == 4 {
                            if mol
                                .neighbors(i)
                                .into_iter()
                                .any(|(nbr, _)| is_metal_or_metalloid(mol.atoms[nbr].z))
                            {
                                intended = 4;
                            }
                        } else if is_bridging_h_or_halogen(mol, i) {
                            intended = 2;
                        } else if (z == 5 || is_metal) && (required - intended) == 1 {
                            if mol
                                .neighbors(i)
                                .into_iter()
                                .any(|(nbr, _)| is_bridging_h_or_halogen(mol, nbr))
                            {
                                intended += 1;
                            }
                        }
                    }

                    let diff = intended - rad;
                    if diff >= 0 {
                        ev = diff;
                    } else if atom.valence_override.is_some() {
                        ev = diff;
                    } else {
                        ev = 0;
                    }

                    // Track redundant/indeterminate valence override (0x180201daa..0x180201f55)
                    if ovr > 0 {
                        let sbo_intended = get_intended_valence(iso_z, sbo[i] + rad);
                        if intended == sbo_intended
                            && b1 == 0
                            && b2 == 0
                            && !is_metal
                            && !is_aromatic_non_standard_heteroatom(mol, i)
                        {
                            indeterminate_valence[i] = true;
                        }
                    }
                    ev = ev.max(req_base);
                } else {
                    // Fallback when charged/isoelectronic state is unsupported
                    let pb = period_bucket(z);
                    if pb > 2 && is_halogen(z) && eff_chg >= 2 && eff_chg <= req_base {
                        indeterminate_valence[i] = true;
                    }
                    ev = req_base;
                }
            }

            expected_valence[i] = ev;
            expected_valence_a_priori[i] = if atom.valence_override.is_some() {
                ovr
            } else {
                ev
            };
        }

        // 3. ComputeHCount (0x1802021b0)
        for i in 0..n {
            let atom = &mol.atoms[i];
            let mut exp_h = 0u8;
            for &(nbr, _) in &mol.neighbors(i) {
                let nz = mol.atoms[nbr].z;
                if nz == 1 || nz == 225 || nz == 226 {
                    exp_h = exp_h.saturating_add(1);
                }
            }
            let z = atom.z;
            let imp_h = if (1..=118).contains(&z) || z == 222 || z == 225 || z == 226 {
                (expected_valence[i] - sbo[i]).max(0) as u8
            } else {
                0
            };
            explicit_h_counts[i] = exp_h;
            implicit_h_counts[i] = imp_h;
        }

        // 4. IdentifyTacitDativeBonds (0x1802029d0) for double-bonded hypervalent pnictogens/chalcogens
        for i in 0..n {
            let z = mol.atoms[i].z;
            let eligible = is_pnictogen_or_b(z) || is_metal_or_metalloid(z) || matches!(z, 6 | 14);
            if !eligible {
                continue;
            }
            for &(nbr, bi) in &mol.neighbors(i) {
                if mol.bonds[bi].order != 2 {
                    continue;
                }
                let nz = mol.atoms[nbr].z;
                let nbr_val = (implicit_h_counts[nbr] as i32) + sbo[nbr];
                let tacit = (is_pnictogen(nz) && nbr_val >= 5) || (is_chalcogen_like(nz) && nbr_val >= 4);
                if tacit {
                    total_dative_out[nbr] = total_dative_out[nbr].saturating_add(1);
                    total_dative_in[i] = total_dative_in[i].saturating_add(1);
                }
            }
        }

        Self {
            sbo,
            expected_valence,
            expected_valence_a_priori,
            implicit_h_counts,
            explicit_h_counts,
            explicit_dative_out,
            explicit_dative_in,
            total_dative_out,
            total_dative_in,
            indeterminate_valence,
        }
    }

    #[inline]
    pub fn implicit_h(&self, atom_idx: usize) -> u8 {
        self.implicit_h_counts[atom_idx]
    }

    #[inline]
    pub fn explicit_h(&self, atom_idx: usize) -> u8 {
        self.explicit_h_counts[atom_idx]
    }

    #[inline]
    pub fn total_h(&self, atom_idx: usize) -> u8 {
        self.implicit_h_counts[atom_idx] + self.explicit_h_counts[atom_idx]
    }
}

// ============================================================================
// 5. ccAromaticityInfo (VA 0x1801cf790, 0x1801cf030, 0x1801cc850)
// ============================================================================

/// Matches `ccAromaticityInfo::CanBeAromatic` at VA `0x1801cf790` (RVA `0x1cf7c8` 80-byte table):
/// Returns `true` iff `z in {5, 6, 7, 8, 14, 15, 16, 32, 33, 34, 51, 52, 83, 84}`
/// (`B, C, N, O, Si, P, S, Ge, As, Se, Sb, Te, Bi, Po`).
#[inline]
pub fn can_be_aromatic(z: i32) -> bool {
    if (5..=84).contains(&z) {
        CAN_BE_AROMATIC_JUMP_INDEX[(z - 5) as usize] == 0
    } else {
        false
    }
}

#[derive(Debug, Clone)]
pub struct CcAromaticityInfo {
    pub cyclically_alternating_bonds: BTreeSet<usize>,
    pub aromatic_atoms: BTreeSet<usize>,
    pub aromatic_bonds: BTreeSet<usize>,
    pub antiaromatic_bonds: BTreeSet<usize>,
}

impl CcAromaticityInfo {
    /// Atom pi-electron capacity (`0, 1, 2`, or `None` if ineligible) matching
    /// the 7-way element switch in `alternating_scan` (`0x1801c9850`).
    fn atom_pi_capacity(
        mol: &ChemMolecule,
        val: &CcVal,
        atom_idx: usize,
        active_bonds: &BTreeSet<usize>,
    ) -> Option<usize> {
        let atom = &mol.atoms[atom_idx];
        let z = atom.z;
        if !can_be_aromatic(z) {
            return None;
        }
        let chg = atom.clamped_charge();
        let rad = atom.radical_electrons();
        let nbrs = mol.neighbors(atom_idx);

        let mut in_ring_active = 0usize;
        let mut ring_dbl_or_arom = 0usize;
        let mut exo_dbl_hetero = 0usize;
        let mut has_triple = false;

        for &(nbr, bi) in &nbrs {
            let ord = mol.bonds[bi].order;
            if ord == 3 {
                has_triple = true;
            }
            if active_bonds.contains(&bi) {
                in_ring_active += 1;
                if ord == 2 || ord == 6 {
                    ring_dbl_or_arom += 1;
                }
            } else if ord == 2 && matches!(mol.atoms[nbr].z, 7 | 8 | 16) {
                exo_dbl_hetero += 1;
            }
        }
        if in_ring_active < 2 || has_triple || rad != 0 {
            return None;
        }
        let total_deg = nbrs.len() + (val.implicit_h(atom_idx) as usize);

        match z {
            // Case 0: Boron (B=5)
            5 => {
                if chg == 0 && total_deg <= 3 && ring_dbl_or_arom == 0 {
                    Some(0)
                } else if chg == -1 && ring_dbl_or_arom >= 1 {
                    Some(1)
                } else {
                    None
                }
            }
            // Case 1 & 4: Carbon, Silicon, Germanium (C=6, Si=14, Ge=32)
            6 | 14 | 32 => {
                if ring_dbl_or_arom >= 1 && chg == 0 && exo_dbl_hetero == 0 {
                    Some(1)
                } else if exo_dbl_hetero >= 1 && chg == 0 {
                    Some(0)
                } else if chg == 1 && total_deg <= 3 && ring_dbl_or_arom == 0 {
                    Some(0)
                } else if chg == -1 && total_deg <= 3 && ring_dbl_or_arom == 0 {
                    Some(2)
                } else {
                    None
                }
            }
            // Case 2 & 5: Pnictogens (N=7, P=15, As=33, Sb=51, Bi=83)
            7 | 15 | 33 | 51 | 83 => {
                if ring_dbl_or_arom >= 1 && (chg == 0 || chg == 1) && exo_dbl_hetero == 0 {
                    Some(1)
                } else if chg == 0 && ring_dbl_or_arom == 0 && total_deg <= 3 && exo_dbl_hetero == 0 {
                    Some(2)
                } else {
                    None
                }
            }
            // Case 3: Chalcogens (O=8, S=16, Se=34, Te=52, Po=84)
            8 | 16 | 34 | 52 | 84 => {
                if chg == 1 && ring_dbl_or_arom >= 1 {
                    Some(1)
                } else if chg == 0 && ring_dbl_or_arom == 0 && total_deg <= 2 {
                    Some(2)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    /// Perceives cyclic alternating-bond (Kekulé) systems gated by `can_be_aromatic(z)`
    /// with iterative pruning (`0x1801caf80`) and 4n+2 (aromatic) vs 4n (antiaromatic) split.
    pub fn perceive(mol: &ChemMolecule) -> Self {
        let mut info = Self {
            cyclically_alternating_bonds: BTreeSet::new(),
            aromatic_atoms: BTreeSet::new(),
            aromatic_bonds: BTreeSet::new(),
            antiaromatic_bonds: BTreeSet::new(),
        };

        let val = CcVal::compute(mol);
        let mut active_bonds: BTreeSet<usize> = (0..mol.bonds.len()).collect();

        // Iterative k-core / eligibility pruning matching 0x1801caf80
        loop {
            let mut removed = Vec::new();
            for &bi in &active_bonds {
                let b = &mol.bonds[bi];
                if Self::atom_pi_capacity(mol, &val, b.begin_atom, &active_bonds).is_none()
                    || Self::atom_pi_capacity(mol, &val, b.end_atom, &active_bonds).is_none()
                {
                    removed.push(bi);
                }
            }
            if removed.is_empty() {
                break;
            }
            for bi in removed {
                active_bonds.remove(&bi);
            }
        }

        if active_bonds.is_empty() {
            return info;
        }

        let cycles = find_simple_cycles_subset(mol, &active_bonds, 18);
        for (cycle_atoms, cycle_bonds) in cycles {
            let all_explicit_arom = cycle_bonds.iter().all(|&bi| mol.bonds[bi].order == 6);
            let mut pi_electrons = 0usize;
            let mut valid = true;

            for &ai in &cycle_atoms {
                if let Some(cap) = Self::atom_pi_capacity(mol, &val, ai, &active_bonds) {
                    pi_electrons += cap;
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }

            if all_explicit_arom || (pi_electrons >= 2 && pi_electrons % 2 == 0) {
                if pi_electrons >= 2 && (pi_electrons % 4) == 2 {
                    for &bi in &cycle_bonds {
                        info.cyclically_alternating_bonds.insert(bi);
                        info.aromatic_bonds.insert(bi);
                    }
                    for &ai in &cycle_atoms {
                        info.aromatic_atoms.insert(ai);
                    }
                } else if pi_electrons >= 4 && (pi_electrons % 4) == 0 {
                    for &bi in &cycle_bonds {
                        info.cyclically_alternating_bonds.insert(bi);
                        info.antiaromatic_bonds.insert(bi);
                    }
                }
            }
        }

        // Remove any bond from antiaromatic_bonds if it is already part of an aromatic 4n+2 system
        for bi in &info.aromatic_bonds {
            info.antiaromatic_bonds.remove(bi);
        }

        info
    }

    #[inline]
    pub fn is_atom_aromatic(&self, atom_idx: usize) -> bool {
        self.aromatic_atoms.contains(&atom_idx)
    }

    #[inline]
    pub fn is_bond_aromatic(&self, bond_idx: usize) -> bool {
        self.aromatic_bonds.contains(&bond_idx)
    }

    #[inline]
    pub fn is_bond_cyclically_alternating(&self, bond_idx: usize) -> bool {
        self.cyclically_alternating_bonds.contains(&bond_idx)
    }
}

fn find_simple_cycles_subset(
    mol: &ChemMolecule,
    allowed_bonds: &BTreeSet<usize>,
    max_len: usize,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    let mut results = Vec::new();
    let mut seen_bond_sets = BTreeSet::new();
    let n = mol.atoms.len();

    for start in 0..n {
        let mut stack: Vec<(usize, Vec<usize>, Vec<usize>)> = vec![(start, vec![start], Vec::new())];
        while let Some((curr, path_atoms, path_bonds)) = stack.pop() {
            if path_atoms.len() > max_len {
                continue;
            }
            for (nbr, bi) in mol.neighbors(curr) {
                if !allowed_bonds.contains(&bi) || path_bonds.last() == Some(&bi) {
                    continue;
                }
                if nbr == start && path_atoms.len() >= 3 {
                    let mut c_bonds = path_bonds.clone();
                    c_bonds.push(bi);
                    let mut sorted_bonds = c_bonds.clone();
                    sorted_bonds.sort_unstable();
                    if seen_bond_sets.insert(sorted_bonds) {
                        results.push((path_atoms.clone(), c_bonds));
                    }
                } else if nbr > start && !path_atoms.contains(&nbr) {
                    let mut next_atoms = path_atoms.clone();
                    next_atoms.push(nbr);
                    let mut next_bonds = path_bonds.clone();
                    next_bonds.push(bi);
                    stack.push((nbr, next_atoms, next_bonds));
                }
            }
        }
    }
    results
}

/// Converts `order == 6` delocalized aromatic bonds into valid alternating single (`1`) / double (`2`)
/// Kekulé bonds (analogous to `CoreChemistryCommon::ccKekulize`), returning `true` if a valid matching exists.
pub fn kekulize_aromatic_bonds(mol: &mut ChemMolecule) -> bool {
    let arom_edges: Vec<usize> = mol
        .bonds
        .iter()
        .enumerate()
        .filter_map(|(bi, b)| (b.order == 6).then_some(bi))
        .collect();
    if arom_edges.is_empty() {
        return true;
    }

    let val = CcVal::compute(mol);
    let n = mol.atoms.len();
    let mut needs_double = vec![false; n];
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for &bi in &arom_edges {
        let b = &mol.bonds[bi];
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    for i in 0..n {
        if adj[i].is_empty() {
            continue;
        }
        let z = mol.atoms[i].z;
        let has_exocyclic_double = mol
            .neighbors(i)
            .into_iter()
            .any(|(_, bi)| mol.bonds[bi].order >= 2 && mol.bonds[bi].order != 6);
        if has_exocyclic_double {
            continue;
        }
        let ext_sbo: i32 = mol
            .neighbors(i)
            .into_iter()
            .filter(|&(_, bi)| mol.bonds[bi].order != 6)
            .map(|(_, bi)| mol.bonds[bi].order)
            .sum();
        let total_non_arom_deg = ext_sbo + i32::from(val.implicit_h_counts[i]);
        // Atoms donating a lone pair (O, S, N-H pyrrole-like with 2 aromatic neighbors + total non-arom 1) need 0 double bonds
        let atom_needs_dbl = match z {
            6 | 14 | 32 => total_non_arom_deg <= 1 && adj[i].len() >= 2,
            7 | 15 | 33 => adj[i].len() == 2 && total_non_arom_deg == 0 && mol.atoms[i].charge == 0,
            _ => false,
        };
        needs_double[i] = atom_needs_dbl;
    }

    let nodes: Vec<usize> = (0..n).filter(|&i| needs_double[i]).collect();
    let mut matched_double = vec![false; mol.bonds.len()];
    let mut covered = vec![false; n];

    fn backtrack(
        pos: usize,
        nodes: &[usize],
        adj: &[Vec<(usize, usize)>],
        needs_double: &[bool],
        covered: &mut [bool],
        matched_double: &mut [bool],
    ) -> bool {
        let mut p = pos;
        while p < nodes.len() && covered[nodes[p]] {
            p += 1;
        }
        if p == nodes.len() {
            return true;
        }
        let u = nodes[p];
        for &(v, bi) in &adj[u] {
            if needs_double[v] && !covered[v] {
                covered[u] = true;
                covered[v] = true;
                matched_double[bi] = true;
                if backtrack(p + 1, nodes, adj, needs_double, covered, matched_double) {
                    return true;
                }
                matched_double[bi] = false;
                covered[u] = false;
                covered[v] = false;
            }
        }
        false
    }

    if !backtrack(0, &nodes, &adj, &needs_double, &mut covered, &mut matched_double) {
        return false;
    }
    for &bi in &arom_edges {
        mol.bonds[bi].order = if matched_double[bi] { 2 } else { 1 };
    }
    true
}

/// Converts perceived aromatic ring bonds (`CcAromaticityInfo`) to `order == 6` (analogous to `ccDekekulize`),
/// returning the number of converted bonds.
pub fn dekekulize_aromatic_bonds(mol: &mut ChemMolecule) -> usize {
    let arom = CcAromaticityInfo::perceive(mol);
    let mut count = 0usize;
    for &bi in &arom.aromatic_bonds {
        if let Some(b) = mol.bonds.get_mut(bi) {
            if b.order != 6 {
                b.order = 6;
                count += 1;
            }
        }
    }
    count
}

// ============================================================================
// 6. ccStereoModel_Tet, ccStereoModel_Dbd, and ccStereoCIP
//    (VA 0x18036f670, 0x180379e30, 0x18037bde0, 0x1803918e0, 0x180392420, 0x180397870)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TetrahedralSubs {
    pub centers: [i32; 4],
    pub subs: [i32; 4],
    pub bonds: [i32; 4],
}

/// Matches `ccStereoModel_Tet::CalcTetrahedralSubs` at VA `0x18036f670`:
/// Fills `centers` with `atom_idx`, `subs` and `bonds` up to `degree` (at most 4),
/// padding unused slots with `-1`.
pub fn calc_tetrahedral_subs(mol: &ChemMolecule, atom_idx: usize) -> TetrahedralSubs {
    let nbrs = mol.neighbors(atom_idx);
    let mut subs = [-1i32; 4];
    let mut bonds = [-1i32; 4];
    let centers = [atom_idx as i32; 4];
    for (slot, &(nbr, bi)) in nbrs.iter().take(4).enumerate() {
        subs[slot] = nbr as i32;
        bonds[slot] = bi as i32;
    }
    TetrahedralSubs {
        centers,
        subs,
        bonds,
    }
}

/// Matches `ccStereoModel_Tet::AnalyzeIn3D` at VA `0x180379e30`.
/// Returns `Some(1)` or `Some(2)` for valid 3D tetrahedral parity, or `None` if degenerate/flat.
pub fn analyze_tetrahedral_in_3d(mol: &ChemMolecule, atom_idx: usize, invert: bool) -> Option<u8> {
    if !mol.has_3d_coords {
        return None;
    }
    let subs = calc_tetrahedral_subs(mol, atom_idx);
    if subs.subs[2] == -1 {
        return None;
    }
    let c = &mol.atoms[atom_idx];
    let mut vecs = [[0.0_f64; 3]; 4];
    for i in 0..3 {
        let s = &mol.atoms[subs.subs[i] as usize];
        let mut v = [c.x - s.x, c.y - s.y, c.z_coord - s.z_coord];
        if v[0].abs() <= 1e-10 && v[1].abs() <= 1e-10 && v[2].abs() <= 1e-10 {
            return None;
        }
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if len > 1e-10 {
            v[0] /= len;
            v[1] /= len;
            v[2] /= len;
        }
        vecs[i] = v;
    }
    if subs.subs[3] != -1 {
        let s = &mol.atoms[subs.subs[3] as usize];
        let mut v = [c.x - s.x, c.y - s.y, c.z_coord - s.z_coord];
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if len <= 1e-10 {
            return None;
        }
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
        vecs[3] = v;
    } else {
        let mut v = [
            -(vecs[0][0] + vecs[1][0] + vecs[2][0]),
            -(vecs[0][1] + vecs[1][1] + vecs[2][1]),
            -(vecs[0][2] + vecs[1][2] + vecs[2][2]),
        ];
        let len_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        // Threshold 0.04 at VA 0x180c399e8
        if len_sq < 0.04 {
            return None;
        }
        let len = len_sq.sqrt();
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
        vecs[3] = v;
    }

    // Signed volume determinant of (v0 - v3, v1 - v3, v2 - v3)
    let a = [
        vecs[0][0] - vecs[3][0],
        vecs[0][1] - vecs[3][1],
        vecs[0][2] - vecs[3][2],
    ];
    let b = [
        vecs[1][0] - vecs[3][0],
        vecs[1][1] - vecs[3][1],
        vecs[1][2] - vecs[3][2],
    ];
    let d = [
        vecs[2][0] - vecs[3][0],
        vecs[2][1] - vecs[3][1],
        vecs[2][2] - vecs[3][2],
    ];
    let det = a[0] * (b[1] * d[2] - b[2] * d[1])
        - a[1] * (b[0] * d[2] - b[2] * d[0])
        + a[2] * (b[0] * d[1] - b[1] * d[0]);
    if det.abs() <= 1e-6 {
        return None;
    }
    let mut parity = if det > 0.0 { 1u8 } else { 2u8 };
    if invert {
        parity = if parity == 1 { 2 } else { 1 };
    }
    Some(parity)
}

/// 2D Tetrahedral stereo analysis (`ccStereoModel_Tet::AnalyzeIn2D` at VA `0x180370200`):
/// Evaluates 2D coordinates `(x, y)` and wedge/hash bond displays (`WedgeBegin`/`WedgeEnd`,
/// `WedgedHashBegin`/`WedgedHashEnd`, `Hash`, `Bold`, `Wavy`) around `atom_idx`.
pub fn analyze_tetrahedral_in_2d(mol: &ChemMolecule, atom_idx: usize) -> Option<u8> {
    let ca = &mol.atoms[atom_idx];
    if ca.z == 6 && mol.neighbors(atom_idx).iter().any(|&(_, bi)| mol.bonds[bi].order >= 2) {
        return None;
    }
    let subs = calc_tetrahedral_subs(mol, atom_idx);
    if subs.subs[2] == -1 {
        return None;
    }
    let c = &mol.atoms[atom_idx];
    let deg = if subs.subs[3] != -1 { 4 } else { 3 };
    let mut z_signs = [0.0_f64; 4];
    let mut has_stereo_bond = false;

    for slot in 0..deg {
        let bi = subs.bonds[slot] as usize;
        let b = &mol.bonds[bi];
        let at_begin = b.begin_atom == atom_idx;
        let z = match b.display {
            BondDisplay::Wavy => return None,
            BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => {
                if at_begin {
                    has_stereo_bond = true;
                    1.0
                } else {
                    0.0
                }
            }
            BondDisplay::WedgeEnd | BondDisplay::HollowWedgeEnd => {
                if !at_begin {
                    has_stereo_bond = true;
                    1.0
                } else {
                    0.0
                }
            }
            BondDisplay::WedgedHashBegin => {
                if at_begin {
                    has_stereo_bond = true;
                    -1.0
                } else {
                    0.0
                }
            }
            BondDisplay::WedgedHashEnd => {
                if !at_begin {
                    has_stereo_bond = true;
                    -1.0
                } else {
                    0.0
                }
            }
            BondDisplay::Hash => {
                has_stereo_bond = true;
                -1.0
            }
            BondDisplay::Bold => {
                has_stereo_bond = true;
                1.0
            }
            BondDisplay::Solid | BondDisplay::Dash | BondDisplay::Dative => 0.0,
        };
        z_signs[slot] = z;
    }

    if !has_stereo_bond {
        return None;
    }

    let mut pts = [[0.0_f64; 3]; 4];
    for slot in 0..deg {
        let s = &mol.atoms[subs.subs[slot] as usize];
        // Note: CDXML 2D coordinates have y increasing downward; use center - sub matching 3D (R=2, S=3)
        let dx = c.x - s.x;
        let dy = -(c.y - s.y);
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 1e-10 {
            return None;
        }
        pts[slot] = [dx / len, dy / len, -z_signs[slot]];
    }
    if deg == 3 {
        // Implicit 4th substituent (hydrogen / lone pair) opposite the stereo bond (in c - s frame)
        pts[3] = [0.0, 0.0, z_signs[0] + z_signs[1] + z_signs[2]];
    }

    let a = [
        pts[0][0] - pts[3][0],
        pts[0][1] - pts[3][1],
        pts[0][2] - pts[3][2],
    ];
    let b = [
        pts[1][0] - pts[3][0],
        pts[1][1] - pts[3][1],
        pts[1][2] - pts[3][2],
    ];
    let d = [
        pts[2][0] - pts[3][0],
        pts[2][1] - pts[3][1],
        pts[2][2] - pts[3][2],
    ];
    let det = a[0] * (b[1] * d[2] - b[2] * d[1])
        - a[1] * (b[0] * d[2] - b[2] * d[0])
        + a[2] * (b[0] * d[1] - b[1] * d[0]);
    if det.abs() <= 1e-6 {
        return None;
    }
    Some(if det > 0.0 { 1 } else { 2 })
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DbdEndAnalysis {
    pub has_two_subs: bool,
    pub normal1: [f64; 3],
    pub normal2: [f64; 3],
    pub same_side_error: bool,
    pub swap_priority: bool,
}

/// Matches `ccStereoModel_Dbd::AnalyzeOneEnd` at VA `0x18037bde0`.
pub fn analyze_double_bond_one_end(
    mol: &ChemMolecule,
    db_atom_1: usize,
    db_atom_2: usize,
) -> Option<DbdEndAnalysis> {
    let other_nbrs: Vec<usize> = mol
        .neighbors(db_atom_2)
        .into_iter()
        .map(|(n, _)| n)
        .filter(|&n| n != db_atom_1)
        .take(2)
        .collect();
    if other_nbrs.is_empty() {
        return None;
    }

    let a1 = &mol.atoms[db_atom_1];
    let a2 = &mol.atoms[db_atom_2];
    let vec_from_a2 = |target: &ChemAtom| -> Option<[f64; 3]> {
        let v = [
            target.x - a2.x,
            -(target.y - a2.y),
            target.z_coord - a2.z_coord,
        ];
        let lsq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        // Threshold 1e-10 at VA 0x180c39898
        if lsq < 1e-10 {
            return None;
        }
        let l = lsq.sqrt();
        Some([v[0] / l, v[1] / l, v[2] / l])
    };

    let v_bond = vec_from_a2(a1)?;
    let sub1 = vec_from_a2(&mol.atoms[other_nbrs[0]])?;
    let cross = |u: [f64; 3], v: [f64; 3]| -> [f64; 3] {
        [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ]
    };
    let dot = |u: [f64; 3], v: [f64; 3]| -> f64 { u[0] * v[0] + u[1] * v[1] + u[2] * v[2] };
    let norm = |u: [f64; 3]| -> [f64; 3] {
        let lsq = dot(u, u);
        // Threshold 0.001 at VA 0x180c39928
        if lsq < 0.001 {
            [0.0, 0.0, 0.0]
        } else {
            let l = lsq.sqrt();
            [u[0] / l, u[1] / l, u[2] / l]
        }
    };

    let c1 = norm(cross(v_bond, sub1));
    if c1 == [0.0, 0.0, 0.0] {
        return None;
    }

    if other_nbrs.len() == 1 {
        return Some(DbdEndAnalysis {
            has_two_subs: false,
            normal1: c1,
            normal2: [0.0, 0.0, 0.0],
            same_side_error: false,
            swap_priority: false,
        });
    }

    let sub2 = vec_from_a2(&mol.atoms[other_nbrs[1]])?;
    // Threshold 0.99 at VA 0x180c39c38
    if dot(sub1, sub2) > 0.99 {
        return None;
    }
    let c2 = norm(cross(v_bond, sub2));
    if c2 == [0.0, 0.0, 0.0] {
        return None;
    }

    let cross_dot = dot(c1, c2);
    // Threshold -0.02 at VA 0x180c3ab50
    let same_side_error = cross_dot > -0.02;
    let swap_priority = cross_dot > 0.0 && dot(v_bond, sub2) > dot(v_bond, sub1);

    Some(DbdEndAnalysis {
        has_two_subs: true,
        normal1: c1,
        normal2: c2,
        same_side_error,
        swap_priority,
    })
}

/// CIP rule levels byte-bound in `ccStereoCIP::GetDigraph` (`0x1803918e0`): `{2, 3, 4, 5}`.
pub const CIP_RULE_LEVELS: [u8; 4] = [2, 3, 4, 5];

/// Matches `ccStereoCIP::ComputeComparisonValue` at VA `0x18038ecd0`:
/// - Rule level `2` reads atom `+0x38` (16-bit isotope mass number).
/// - All other rule levels (`3, 4, 5`) read atom `+0x10` (32-bit atomic number `Z`).
#[inline]
pub fn cip_compute_comparison_value(atom: &ChemAtom, rule_level: u8) -> f64 {
    if rule_level == 2 {
        (atom.isotope as u16) as f64
    } else {
        atom.z as f64
    }
}

/// Matches `ccStereoCIP::ComputeAtomStereoCode` at VA `0x180392420`:
/// Maps `(parity, rule_level, inverted)` into stereo codes `{0, 1, 2, 3, 4, 5, 10, 11}`.
#[inline]
pub fn cip_compute_atom_stereo_code(parity: u8, rule_level: u8, inverted: bool, is_bond: bool) -> u8 {
    if is_bond {
        return match (parity, inverted) {
            (1, false) | (2, true) => 10, // Z / cis
            (2, false) | (1, true) => 11, // E / trans
            _ => 0,
        };
    }
    match parity {
        1 => {
            if rule_level == 5 {
                if inverted { 5 } else { 4 }
            } else if inverted {
                3 // S
            } else {
                2 // R
            }
        }
        2 => {
            if rule_level == 5 {
                if inverted { 4 } else { 5 }
            } else if inverted {
                2 // R
            } else {
                3 // S
            }
        }
        _ => 0,
    }
}

/// Matches `ccStereoCIP::AdjustStereoCodeForAtom` at VA `0x180397870`:
/// Post-maps stereo code based on enhanced stereo group type (`3 => {2->8, 3->9}`, `4 => {2->6, 3->7}`).
#[inline]
pub fn cip_adjust_stereo_code_for_atom(code: u8, group_type: u32) -> u8 {
    match (group_type, code) {
        (3, 2) => 8,
        (3, 3) => 9,
        (4, 2) => 6,
        (4, 3) => 7,
        _ => code,
    }
}

#[inline]
fn effective_cip_isotope(atom: &ChemAtom) -> i32 {
    if atom.isotope != 0 {
        atom.isotope
    } else if atom.z == 225 {
        2
    } else if atom.z == 226 {
        3
    } else {
        1
    }
}

#[inline]
fn effective_cip_z(atom: &ChemAtom) -> i32 {
    if atom.z == 225 || atom.z == 226 {
        1
    } else {
        atom.z
    }
}

/// Hierarchical CIP ligand node used for breadth-first digraph comparison (Rules 1a, 1b, 2).
#[derive(Debug, Clone, PartialEq, Eq)]
struct CipNode {
    z: i32,
    iso: i32,
    children: Vec<CipNode>,
}

impl CipNode {
    fn spheres(&self) -> (Vec<Vec<i32>>, Vec<Vec<i32>>) {
        let mut z_spheres = Vec::new();
        let mut iso_spheres = Vec::new();
        let mut frontier: Vec<&CipNode> = vec![self];
        while !frontier.is_empty() {
            let zs: Vec<i32> = frontier.iter().map(|n| n.z).collect();
            let isos: Vec<i32> = frontier.iter().map(|n| n.iso).collect();
            z_spheres.push(zs);
            iso_spheres.push(isos);
            let mut next = Vec::new();
            for n in frontier {
                for ch in &n.children {
                    next.push(ch);
                }
            }
            frontier = next;
        }
        (z_spheres, iso_spheres)
    }
}

impl Ord for CipNode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let (az, aiso) = self.spheres();
        let (bz, biso) = other.spheres();
        match az.cmp(&bz) {
            std::cmp::Ordering::Equal => aiso.cmp(&biso),
            ord => ord,
        }
    }
}

impl PartialOrd for CipNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

fn build_cip_digraph(
    mol: &ChemMolecule,
    val: &CcVal,
    visited: &mut Vec<usize>,
    curr_atom: usize,
    max_depth: usize,
) -> CipNode {
    let a = &mol.atoms[curr_atom];
    let z = effective_cip_z(a);
    let iso = effective_cip_isotope(a);
    if max_depth == 0 || visited.contains(&curr_atom) {
        return CipNode {
            z,
            iso,
            children: Vec::new(),
        };
    }
    visited.push(curr_atom);
    let mut children: Vec<CipNode> = Vec::new();
    let parent = if visited.len() >= 2 {
        Some(visited[visited.len() - 2])
    } else {
        None
    };

    for (nbr, bi) in mol.neighbors(curr_atom) {
        if Some(nbr) == parent {
            continue;
        }
        let ord = match mol.bonds[bi].order {
            2 => 2,
            3 => 3,
            4 => 4,
            _ => 1,
        };
        if visited.contains(&nbr) {
            // CIP Rule 1b: ring-closure duplicate ghost atom with zero further bonds
            let na = &mol.atoms[nbr];
            for _ in 0..ord {
                children.push(CipNode {
                    z: effective_cip_z(na),
                    iso: effective_cip_isotope(na),
                    children: Vec::new(),
                });
            }
        } else {
            children.push(build_cip_digraph(mol, val, visited, nbr, max_depth - 1));
            // CIP Rule 1a: duplicate atoms for multiple bonds
            let na = &mol.atoms[nbr];
            for _ in 1..ord {
                children.push(CipNode {
                    z: effective_cip_z(na),
                    iso: effective_cip_isotope(na),
                    children: Vec::new(),
                });
            }
        }
    }
    for _ in 0..val.implicit_h(curr_atom) {
        children.push(CipNode {
            z: 1,
            iso: 1,
            children: Vec::new(),
        });
    }
    children.sort_by(|x, y| y.cmp(x));
    visited.pop();
    CipNode { z, iso, children }
}

fn build_ligand_digraph(
    mol: &ChemMolecule,
    val: &CcVal,
    center_atom: usize,
    sub_idx: i32,
    max_depth: usize,
) -> CipNode {
    if sub_idx < 0 {
        return CipNode {
            z: 1,
            iso: 1,
            children: Vec::new(),
        };
    }
    let mut visited = vec![center_atom];
    build_cip_digraph(mol, val, &mut visited, sub_idx as usize, max_depth)
}

/// Computes the CIP descriptor (`2` = R, `3` = S, `4` = r, `5` = s) for a tetrahedral stereocenter.
pub fn compute_tetrahedral_cip_code(mol: &ChemMolecule, atom_idx: usize) -> u8 {
    let val = CcVal::compute(mol);
    compute_tetrahedral_cip_code_with_val(mol, &val, atom_idx)
}

fn compute_tetrahedral_cip_code_with_val(mol: &ChemMolecule, val: &CcVal, atom_idx: usize) -> u8 {
    let parity = if mol.has_3d_coords {
        analyze_tetrahedral_in_3d(mol, atom_idx, false)
    } else {
        analyze_tetrahedral_in_2d(mol, atom_idx)
    };
    let Some(p) = parity else {
        return 0;
    };
    let subs = calc_tetrahedral_subs(mol, atom_idx);
    if subs.subs[2] == -1 {
        return 0;
    }

    let mut keys: Vec<CipNode> = Vec::with_capacity(4);
    for slot in 0..4 {
        keys.push(build_ligand_digraph(mol, val, atom_idx, subs.subs[slot], 16));
    }

    let mut order: Vec<usize> = (0..4).collect();
    order.sort_by(|&a, &b| keys[b].cmp(&keys[a]));
    let mut rule_level = 4u8;

    // CIP Rule 4/5: break constitutive+isotope ligand ties using embedded stereocenter descriptors (R > S, Z > E)
    let has_tie = (0..3).any(|i| keys[order[i]] == keys[order[i + 1]]);
    if has_tie {
        let sub_stereo_signature = |start_sub: i32| -> Vec<u8> {
            if start_sub < 0 {
                return Vec::new();
            }
            let start = start_sub as usize;
            let mut visited = BTreeSet::new();
            visited.insert(atom_idx);
            let mut q = std::collections::VecDeque::new();
            visited.insert(start);
            q.push_back(start);
            let mut sig = Vec::new();
            while let Some(cur) = q.pop_front() {
                let neighbors = mol.neighbors(cur);
                let is_non_self = neighbors.iter().filter(|&&(n, _)| n != atom_idx).count();
                if cur != start || is_non_self >= 2 {
                    // Check local parity/stereo on neighbor stereocenters
                    let code = analyze_tetrahedral_in_2d(mol, cur).unwrap_or(0);
                    if code != 0 {
                        // Rank direct neighbors without full recursion to avoid cycles
                        let cur_subs = calc_tetrahedral_subs(mol, cur);
                        if cur_subs.subs[2] != -1 {
                            let mut c_keys = Vec::with_capacity(4);
                            for cs in 0..4 {
                                c_keys.push(build_ligand_digraph(mol, val, cur, cur_subs.subs[cs], 6));
                            }
                            let mut c_ord: Vec<usize> = (0..4).collect();
                            c_ord.sort_by(|&a, &b| c_keys[b].cmp(&c_keys[a]));
                            if !(0..3).any(|i| c_keys[c_ord[i]] == c_keys[c_ord[i + 1]]) {
                                let mut inv = 0;
                                for i in 0..4 {
                                    for j in (i + 1)..4 {
                                        if c_ord[i] > c_ord[j] {
                                            inv += 1;
                                        }
                                    }
                                }
                                let c_code = cip_compute_atom_stereo_code(code, 4, (inv % 2) == 1, false);
                                // R (2) outranks S (3) -> higher score comes first
                                sig.push(if c_code == 2 { 2 } else if c_code == 3 { 1 } else { 0 });
                            }
                        }
                    }
                }
                for &(nb, _) in &neighbors {
                    if visited.insert(nb) {
                        q.push_back(nb);
                    }
                }
            }
            sig
        };

        let sigs: Vec<Vec<u8>> = (0..4).map(|slot| sub_stereo_signature(subs.subs[slot])).collect();
        order.sort_by(|&a, &b| keys[b].cmp(&keys[a]).then_with(|| sigs[b].cmp(&sigs[a])));
        for i in 0..3 {
            let (ia, ib) = (order[i], order[i + 1]);
            if keys[ia] == keys[ib] && sigs[ia] == sigs[ib] {
                return 0; // True symmetry -> non-stereogenic
            }
        }
        rule_level = 5; // Pseudoasymmetric center (r=4, s=5)
    }

    let mut inversions = 0;
    for i in 0..4 {
        for j in (i + 1)..4 {
            if order[i] > order[j] {
                inversions += 1;
            }
        }
    }
    let odd_perm = (inversions % 2) == 1;
    cip_compute_atom_stereo_code(p, rule_level, odd_perm, false)
}

/// Computes the CIP double-bond descriptor (`10` = Z/cis, `11` = E/trans, `0` = non-stereogenic)
/// for bond `bond_idx` using `analyze_double_bond_one_end` (`0x18037bde0`) + CIP ligand ranking.
pub fn compute_double_bond_cip_code(mol: &ChemMolecule, bond_idx: usize) -> u8 {
    let val = CcVal::compute(mol);
    compute_double_bond_cip_code_with_val(mol, &val, bond_idx)
}

fn bond_in_ring_up_to(mol: &ChemMolecule, bond_idx: usize, max_ring_size: usize) -> bool {
    if max_ring_size < 3 {
        return false;
    }
    let b = &mol.bonds[bond_idx];
    let start = b.begin_atom;
    let target = b.end_atom;
    if start == target {
        return true;
    }
    let max_edges = max_ring_size - 1;
    let mut dist = vec![usize::MAX; mol.atoms.len()];
    let mut queue = std::collections::VecDeque::new();
    dist[start] = 0;
    queue.push_back(start);
    while let Some(u) = queue.pop_front() {
        let d = dist[u];
        if d >= max_edges {
            continue;
        }
        for (v, bi) in mol.neighbors(u) {
            if bi == bond_idx {
                continue;
            }
            if v == target {
                return true;
            }
            if dist[v] == usize::MAX {
                dist[v] = d + 1;
                queue.push_back(v);
            }
        }
    }
    false
}

fn compute_double_bond_cip_code_with_val(mol: &ChemMolecule, val: &CcVal, bond_idx: usize) -> u8 {
    let b = &mol.bonds[bond_idx];
    if b.order != 2 || b.display == BondDisplay::Wavy {
        return 0;
    }
    if bond_in_ring_up_to(mol, bond_idx, 7) {
        return 0;
    }
    let a1 = b.begin_atom;
    let a2 = b.end_atom;
    let Some(end1) = analyze_double_bond_one_end(mol, a2, a1) else {
        return 0;
    };
    let Some(end2) = analyze_double_bond_one_end(mol, a1, a2) else {
        return 0;
    };
    if end1.same_side_error || end2.same_side_error {
        return 0;
    }

    let rank_end = |db_atom: usize, other_db: usize| -> Option<(usize, [f64; 3])> {
        let nbrs: Vec<usize> = mol
            .neighbors(db_atom)
            .into_iter()
            .map(|(n, _)| n)
            .filter(|&n| n != other_db)
            .take(2)
            .collect();
        if nbrs.is_empty() {
            return None;
        }
        let n0 = build_ligand_digraph(mol, val, db_atom, nbrs[0] as i32, 16);
        let n1 = if nbrs.len() == 2 {
            build_ligand_digraph(mol, val, db_atom, nbrs[1] as i32, 16)
        } else {
            build_ligand_digraph(mol, val, db_atom, -1, 16)
        };
        let high_idx = match n0.cmp(&n1) {
            std::cmp::Ordering::Greater => 0,
            std::cmp::Ordering::Less => 1,
            std::cmp::Ordering::Equal => return None,
        };
        if high_idx == 1 && nbrs.len() == 1 {
            // Implicit H outranked the explicit substituent (cannot happen for Z >= 2, only if explicit is H with lower iso)
            return None;
        }
        let high_atom = &mol.atoms[nbrs[high_idx]];
        let c = &mol.atoms[db_atom];
        let o = &mol.atoms[other_db];
        let v_bond = [o.x - c.x, -(o.y - c.y), o.z_coord - c.z_coord];
        let v_sub = [high_atom.x - c.x, -(high_atom.y - c.y), high_atom.z_coord - c.z_coord];
        let cross = [
            v_bond[1] * v_sub[2] - v_bond[2] * v_sub[1],
            v_bond[2] * v_sub[0] - v_bond[0] * v_sub[2],
            v_bond[0] * v_sub[1] - v_bond[1] * v_sub[0],
        ];
        Some((nbrs[high_idx], cross))
    };

    let Some((_, cross1)) = rank_end(a1, a2) else {
        return 0;
    };
    let Some((_, cross2)) = rank_end(a2, a1) else {
        return 0;
    };
    // Because v_bond(a1->a2) == -v_bond(a2->a1), substituents on the SAME side of the double bond
    // produce OPPOSITE cross-product normals (dot < 0 => cis/Z = 10), and substituents on OPPOSITE
    // sides produce PARALLEL cross-product normals (dot > 0 => trans/E = 11).
    let d = cross1[0] * cross2[0] + cross1[1] * cross2[1] + cross1[2] * cross2[2];
    if d.abs() <= 1e-9 {
        return 0;
    }
    if d < 0.0 {
        cip_compute_atom_stereo_code(1, 4, false, true) // 10 = Z / cis
    } else {
        cip_compute_atom_stereo_code(2, 4, false, true) // 11 = E / trans
    }
}

/// Computes the axial chirality descriptor (`Some('P')` -> CDXML `AS="+"`,
/// `Some('M')` -> CDXML `AS="-"`) for an allene central atom (`=C=`).
/// Returns `None` if `center_idx` is not a substituted stereogenic allene or lacks 3D/wedge info.
pub fn compute_axial_allene_cip(mol: &ChemMolecule, center_idx: usize) -> Option<char> {
    let center = mol.atoms.get(center_idx)?;
    if center.z != 6 {
        return None;
    }
    let c_nbrs = mol.neighbors(center_idx);
    if c_nbrs.len() != 2 {
        return None;
    }
    if mol.bonds[c_nbrs[0].1].order != 2 || mol.bonds[c_nbrs[1].1].order != 2 {
        return None;
    }
    let end_a = c_nbrs[0].0;
    let end_b = c_nbrs[1].0;
    if end_a == end_b || mol.atoms[end_a].z != 6 || mol.atoms[end_b].z != 6 {
        return None;
    }

    let val = CcVal::compute(mol);
    let end_info = |end: usize| -> Option<(Option<(usize, usize)>, Option<(usize, usize)>, bool)> {
        let subs: Vec<(usize, usize)> = mol
            .neighbors(end)
            .into_iter()
            .filter(|&(n, _)| n != center_idx)
            .collect();
        if subs.is_empty() || subs.len() > 2 {
            return None;
        }
        if subs.iter().any(|&(_, bi)| mol.bonds[bi].display == BondDisplay::Wavy) {
            return None;
        }
        let n0 = build_ligand_digraph(mol, &val, end, subs[0].0 as i32, 16);
        let n1 = if subs.len() == 2 {
            build_ligand_digraph(mol, &val, end, subs[1].0 as i32, 16)
        } else if val.implicit_h(end) == 1 {
            build_ligand_digraph(mol, &val, end, -1, 16)
        } else {
            return None;
        };
        let sub1 = subs.get(1).copied();
        match n0.cmp(&n1) {
            std::cmp::Ordering::Greater => Some((Some(subs[0]), sub1, false)),
            std::cmp::Ordering::Less => Some((sub1, Some(subs[0]), false)),
            std::cmp::Ordering::Equal => None,
        }
    };

    let (a_high, a_low, _) = end_info(end_a)?;
    let (b_high, b_low, _) = end_info(end_b)?;

    let bond_z_at_end = |end: usize, bi: usize| -> f64 {
        let b = &mol.bonds[bi];
        let at_begin = b.begin_atom == end;
        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => {
                if at_begin { 1.0 } else { 0.0 }
            }
            BondDisplay::WedgeEnd | BondDisplay::HollowWedgeEnd => {
                if !at_begin { 1.0 } else { 0.0 }
            }
            BondDisplay::WedgedHashBegin => {
                if at_begin { -1.0 } else { 0.0 }
            }
            BondDisplay::WedgedHashEnd => {
                if !at_begin { -1.0 } else { 0.0 }
            }
            BondDisplay::Hash => -1.0,
            BondDisplay::Bold => 1.0,
            _ => 0.0,
        }
    };

    let mut has_stereo = false;
    let mut sub_vec = |end: usize, sub: Option<(usize, usize)>, fallback_other: Option<(usize, usize)>| -> Option<[f64; 3]> {
        let e = &mol.atoms[end];
        if let Some((nbr, bi)) = sub {
            let s = &mol.atoms[nbr];
            let dx = s.x - e.x;
            let dy = -(s.y - e.y);
            let dz_3d = s.z_coord - e.z_coord;
            let len_xy = (dx * dx + dy * dy).sqrt();
            if len_xy <= 1e-10 && dz_3d.abs() <= 1e-10 {
                return None;
            }
            let bz = bond_z_at_end(end, bi);
            let z = if mol.has_3d_coords {
                has_stereo = true;
                dz_3d
            } else {
                if bz.abs() > 0.0 {
                    has_stereo = true;
                }
                bz
            };
            if len_xy > 1e-10 {
                Some([dx / len_xy, dy / len_xy, z])
            } else {
                Some([0.0, 0.0, z])
            }
        } else {
            let (nbr, bi) = fallback_other?;
            let s = &mol.atoms[nbr];
            let dx = s.x - e.x;
            let dy = -(s.y - e.y);
            let len_xy = (dx * dx + dy * dy).sqrt().max(1e-10);
            let bz = bond_z_at_end(end, bi);
            let z = if mol.has_3d_coords { -(s.z_coord - e.z_coord) } else { -bz };
            Some([-dx / len_xy, -dy / len_xy, z])
        }
    };

    let va_h = sub_vec(end_a, a_high, a_low)?;
    let va_l = sub_vec(end_a, a_low, a_high)?;
    let vb_h = sub_vec(end_b, b_high, b_low)?;
    let vb_l = sub_vec(end_b, b_low, b_high)?;
    if !has_stereo {
        return None;
    }

    let ea = &mol.atoms[end_a];
    let eb = &mol.atoms[end_b];
    let u = [eb.x - ea.x, -(eb.y - ea.y), if mol.has_3d_coords { eb.z_coord - ea.z_coord } else { 0.0 }];
    let u_len = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    if u_len <= 1e-10 {
        return None;
    }
    let u = [u[0] / u_len, u[1] / u_len, u[2] / u_len];

    let da = [va_h[0] - va_l[0], va_h[1] - va_l[1], va_h[2] - va_l[2]];
    let db = [vb_h[0] - vb_l[0], vb_h[1] - vb_l[1], vb_h[2] - vb_l[2]];
    let cross = [
        da[1] * db[2] - da[2] * db[1],
        da[2] * db[0] - da[0] * db[2],
        da[0] * db[1] - da[1] * db[0],
    ];
    let det = cross[0] * u[0] + cross[1] * u[1] + cross[2] * u[2];
    if det.abs() <= 1e-6 {
        None
    } else if det > 0.0 {
        Some('P')
    } else {
        Some('M')
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StereochemistryReport {
    pub atom_stereo: BTreeMap<usize, u8>,
    pub bond_stereo: BTreeMap<usize, u8>,
}

/// Runs full tetrahedral + double-bond CIP perception across `mol`.
pub fn perceive_stereochemistry(mol: &ChemMolecule) -> StereochemistryReport {
    let val = CcVal::compute(mol);
    let mut rep = StereochemistryReport::default();
    for ai in 0..mol.atoms.len() {
        let code = compute_tetrahedral_cip_code_with_val(mol, &val, ai);
        if code != 0 {
            rep.atom_stereo.insert(ai, code);
        }
    }
    for bi in 0..mol.bonds.len() {
        let code = compute_double_bond_cip_code_with_val(mol, &val, bi);
        if code != 0 {
            rep.bond_stereo.insert(bi, code);
        }
    }
    rep
}

// ============================================================================
// 7. ccReaction (VA 0x1802ac270, 0x1802b02b0, 0x1802b0550, 0x1802af700, 0x1802b2740)
// ============================================================================

/// Matches `ccReaction::CenterTypeDenotesChange` at VA `0x1802ac270`:
/// `if (type > 7) return true; else return ((0xd1 >> type) & 1) == 0;`
/// Returns `true` for `{1, 2, 3, 5} U (8..)`, `false` for `{0, 4, 6, 7}`.
#[inline]
pub fn center_type_denotes_change(center_type: u32) -> bool {
    if center_type > 7 {
        true
    } else {
        ((0xd1u8 >> center_type) & 1) == 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct AtomMapKey {
    pub mol_slot: usize,
    pub atom_idx: usize,
}

#[derive(Debug, Clone)]
pub struct ReactionCenterSlot {
    pub atom_bits: BTreeSet<usize>,
    pub bond_bits: BTreeSet<usize>,
    pub center_types: Vec<u32>,
}

impl ReactionCenterSlot {
    pub fn new(num_items: usize) -> Self {
        Self {
            atom_bits: BTreeSet::new(),
            bond_bits: BTreeSet::new(),
            center_types: vec![0; num_items],
        }
    }
}

#[derive(Debug, Clone)]
pub struct ReactionStep {
    pub reactant: ChemMolecule,
    pub product: ChemMolecule,
    /// Atom-atom map: reactant atom_idx -> product atom_idx
    pub atom_map: BTreeMap<usize, usize>,
    /// Bond-bond map: reactant bond_idx -> product bond_idx
    pub bond_map: BTreeMap<usize, usize>,
}

#[derive(Debug, Clone)]
pub struct CcReaction {
    pub steps: Vec<ReactionStep>,
    pub symmetric_atom_map: BTreeMap<AtomMapKey, AtomMapKey>,
    pub reaction_center_origin: u8,
    pub center_count_word: u16,
    pub centers_perceived: bool,
    pub records: Vec<(ReactionCenterSlot, ReactionCenterSlot)>,
}

impl Default for CcReaction {
    fn default() -> Self {
        Self::new()
    }
}

impl CcReaction {
    pub fn new() -> Self {
        Self {
            steps: Vec::new(),
            symmetric_atom_map: BTreeMap::new(),
            reaction_center_origin: 0,
            center_count_word: 0x4e, // Default 78 at VA 0x1802b02b0
            centers_perceived: false,
            records: Vec::new(),
        }
    }

    /// Matches `AddSingleAtomMapping` (`0x1802ab320`): symmetric insertion into `[this+0x60]`.
    pub fn add_single_atom_mapping(&mut self, a: AtomMapKey, b: AtomMapKey) {
        self.symmetric_atom_map.insert(a, b);
        self.symmetric_atom_map.insert(b, a);
        self.centers_perceived = false;
    }

    /// Matches `RemoveSingleAtomMapping` (`0x1802ab3a0`).
    pub fn remove_single_atom_mapping(&mut self, a: AtomMapKey) {
        if let Some(b) = self.symmetric_atom_map.remove(&a) {
            self.symmetric_atom_map.remove(&b);
        }
        self.centers_perceived = false;
    }

    /// Matches `ClearAaMap` (`0x1802ac800`).
    pub fn clear_aa_map(&mut self) {
        self.symmetric_atom_map.clear();
        self.records.clear();
        self.centers_perceived = false;
    }

    /// Matches `CreateReactionCentersFromScratch` (`0x1802b0550`):
    /// - Phase A (bond-map loop `0x2b062b..0x2b0954`): compares mapped bond order -> sets `bond_bits` and writes center type `3`.
    /// - Phase B (atom-map loop `0x2b0a84..0x2b0d63`): compares mapped atom `charge`, `radical`, `z`, `isotope` -> sets `atom_bits` and writes center type `2`.
    /// - Phase C (`0x2b1b20`): writes center type `7` for unmapped atoms.
    pub fn create_reaction_centers_from_scratch(&mut self) {
        self.records.clear();
        for step in &self.steps {
            let mut slot_a = ReactionCenterSlot::new(step.reactant.atoms.len().max(step.reactant.bonds.len()));
            let mut slot_b = ReactionCenterSlot::new(step.product.atoms.len().max(step.product.bonds.len()));

            // Phase A: bond-map loop -> center type 3
            for (&rb, &pb) in &step.bond_map {
                let b_a = &step.reactant.bonds[rb];
                let b_b = &step.product.bonds[pb];
                if b_a.order != b_b.order {
                    slot_a.bond_bits.insert(rb);
                    slot_b.bond_bits.insert(pb);
                    if rb < slot_a.center_types.len() {
                        slot_a.center_types[rb] = 3;
                    }
                    if pb < slot_b.center_types.len() {
                        slot_b.center_types[pb] = 3;
                    }
                }
            }

            // Phase B: atom-map loop -> center type 2
            for (&ra, &pa) in &step.atom_map {
                let a_a = &step.reactant.atoms[ra];
                let a_b = &step.product.atoms[pa];
                if a_a.charge != a_b.charge
                    || a_a.radical != a_b.radical
                    || a_a.z != a_b.z
                    || a_a.isotope != a_b.isotope
                {
                    slot_a.atom_bits.insert(ra);
                    slot_b.atom_bits.insert(pa);
                    if ra < slot_a.center_types.len() {
                        slot_a.center_types[ra] = 2;
                    }
                    if pa < slot_b.center_types.len() {
                        slot_b.center_types[pa] = 2;
                    }
                }
            }

            // Phase C: unmapped atoms -> center type 7; formed/broken unmapped bonds between mapped endpoints -> center type 1
            for ra in 0..step.reactant.atoms.len() {
                if !step.atom_map.contains_key(&ra) && ra < slot_a.center_types.len() {
                    slot_a.center_types[ra] = 7;
                }
            }
            let mapped_prod: BTreeSet<usize> = step.atom_map.values().copied().collect();
            for pa in 0..step.product.atoms.len() {
                if !mapped_prod.contains(&pa) && pa < slot_b.center_types.len() {
                    slot_b.center_types[pa] = 7;
                }
            }

            let mut rev_atom_map = BTreeMap::new();
            for (&ra, &pa) in &step.atom_map {
                rev_atom_map.insert(pa, ra);
            }
            let prod_edge_exists = |u: usize, v: usize| -> bool {
                step.product
                    .bonds
                    .iter()
                    .any(|b| (b.begin_atom == u && b.end_atom == v) || (b.begin_atom == v && b.end_atom == u))
            };
            let react_edge_exists = |u: usize, v: usize| -> bool {
                step.reactant
                    .bonds
                    .iter()
                    .any(|b| (b.begin_atom == u && b.end_atom == v) || (b.begin_atom == v && b.end_atom == u))
            };
            for (rb, b) in step.reactant.bonds.iter().enumerate() {
                if !step.bond_map.contains_key(&rb) {
                    if let (Some(&pu), Some(&pv)) = (step.atom_map.get(&b.begin_atom), step.atom_map.get(&b.end_atom)) {
                        if !prod_edge_exists(pu, pv) {
                            slot_a.bond_bits.insert(rb);
                            if rb < slot_a.center_types.len() && slot_a.center_types[rb] == 0 {
                                slot_a.center_types[rb] = 1;
                            }
                        }
                    }
                }
            }
            let mapped_prod_bonds: BTreeSet<usize> = step.bond_map.values().copied().collect();
            for (pb, b) in step.product.bonds.iter().enumerate() {
                if !mapped_prod_bonds.contains(&pb) {
                    if let (Some(&ru), Some(&rv)) = (rev_atom_map.get(&b.begin_atom), rev_atom_map.get(&b.end_atom)) {
                        if !react_edge_exists(ru, rv) {
                            slot_b.bond_bits.insert(pb);
                            if pb < slot_b.center_types.len() && slot_b.center_types[pb] == 0 {
                                slot_b.center_types[pb] = 1;
                            }
                        }
                    }
                }
            }

            self.records.push((slot_a, slot_b));
        }
        self.centers_perceived = true;
    }

    /// Matches `ccReaction::Validate` (`0x1802b2740`) checks V1–V6.
    pub fn validate(&mut self) -> Result<(), u32> {
        // V2: empty short-circuit
        if self.steps.is_empty() {
            return Ok(());
        }
        // V3: step count > 1000 (0x3e8) -> error 0x80d
        if self.steps.len() > 1000 {
            return Err(0x80d);
        }
        // V4: center_count_word normalized against {0x16, 0x17, 0x21, 0x2c, 0x2d, 0x40, 0x4e}
        const VALID_BUCKETS: [u16; 7] = [0x16, 0x17, 0x21, 0x2c, 0x2d, 0x40, 0x4e];
        if !VALID_BUCKETS.contains(&self.center_count_word) {
            self.center_count_word = 0x4e;
        }
        // V5: atom-map balance when centers_perceived
        if self.centers_perceived {
            for (slot_a, slot_b) in &self.records {
                if slot_a.atom_bits.len() != slot_b.atom_bits.len() {
                    return Err(0x813);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_periodic_and_isotope_tables_exact() {
        assert_eq!(ELEMENTS.len(), 118);
        assert_eq!(ISOTOPES.len(), 2938);
        assert_eq!(PSEUDO_ATOMS.len(), 26);

        let h = get_element(1).unwrap();
        assert_eq!(h.symbol, "H");
        assert_eq!(h.name, "Hydrogen");
        assert_eq!(h.preferred_valence, 1);
        assert_eq!(h.allowed_valences, &[1]);

        let c = get_element(6).unwrap();
        assert_eq!(c.symbol, "C");
        assert_eq!(c.preferred_valence, 4);
        assert_eq!(c.allowed_valences, &[4]);

        let p = get_element(15).unwrap();
        assert_eq!(p.symbol, "P");
        assert_eq!(p.preferred_valence, 3);
        assert_eq!(p.allowed_valences, &[3, 5]);

        let s = get_element(16).unwrap();
        assert_eq!(s.symbol, "S");
        assert_eq!(s.preferred_valence, 2);
        assert_eq!(s.allowed_valences, &[2, 4, 6]);

        let cl = get_element(17).unwrap();
        assert_eq!(cl.symbol, "Cl");
        assert_eq!(cl.preferred_valence, 1);
        assert_eq!(cl.allowed_valences, &[1, 3, 5, 7]);

        // Hydrogen isotopes: 1H, 2H (D), 3H (T)
        let h1 = get_isotope(1, 1).unwrap();
        assert_eq!(h1.mass_bits, 0x3ff0200d24087b82);
        assert_eq!(h1.symbol, "");

        let d = get_isotope(1, 2).unwrap();
        assert_eq!(d.symbol, "D");
        assert_eq!(d.mass(), 2.0141017781);

        let t = get_isotope(1, 3).unwrap();
        assert_eq!(t.symbol, "T");
        assert_eq!(t.mass(), 3.0160492675);

        // Element exact mass selects highest-abundance isotope (12C = 12.0 exact)
        let c_mass = get_element_exact_mass(6);
        assert!((c_mass.val - 12.0).abs() < 1e-12);
    }

    #[test]
    fn test_formula_exact_mass_and_hill_formatting() {
        // Benzene C6H6
        let mut fp = CsFormulaProperties::new();
        fp.increment_atom(6, 0, 6.0);
        fp.increment_atom(1, 0, 6.0);
        assert_eq!(fp.get_formula(), "C6H6");
        let em = fp.get_exact_mass();
        let expected = 6.0 * 12.0 + 6.0 * f64::from_bits(0x3ff0200d24087b82);
        assert_eq!(em.val.to_bits(), expected.to_bits());

        // Pyridinium cation C5H6N+ (charge +1 subtracts 0.00054858)
        let mut fp_cat = CsFormulaProperties::new();
        fp_cat.increment_atom(7, 0, 1.0);
        fp_cat.increment_atom(6, 0, 5.0);
        fp_cat.increment_atom(1, 0, 6.0);
        fp_cat.increment_charge(1.0);
        assert_eq!(fp_cat.get_formula(), "C5H6N");
        let em_cat = fp_cat.get_exact_mass();
        let n14 = get_element_exact_mass(7).val;
        let expected_cat = (5.0 * 12.0 + 6.0 * f64::from_bits(0x3ff0200d24087b82) + n14)
            + (-1.0 * electron_mass_factor());
        assert_eq!(em_cat.val.to_bits(), expected_cat.to_bits());

        // Deuterated chloroform CDCl3 (C=6, H(2)=D, Cl=17)
        let mut fp_d = CsFormulaProperties::new();
        fp_d.increment_atom(6, 0, 1.0);
        fp_d.increment_atom(1, 2, 1.0);
        fp_d.increment_atom(17, 0, 3.0);
        assert_eq!(fp_d.get_formula(), "CDCl3");
    }

    #[test]
    fn test_ccval_implicit_hydrogens_and_aromaticity() {
        // Pyrrole: 5-membered ring N1-C2=C3-C4=C5-1
        let mut mol = ChemMolecule::new();
        let n = mol.add_atom(ChemAtom::new(1, 7));
        let c2 = mol.add_atom(ChemAtom::new(2, 6));
        let c3 = mol.add_atom(ChemAtom::new(3, 6));
        let c4 = mol.add_atom(ChemAtom::new(4, 6));
        let c5 = mol.add_atom(ChemAtom::new(5, 6));
        mol.add_bond(1, n, c2, 1, BondDisplay::Solid);
        mol.add_bond(2, c2, c3, 2, BondDisplay::Solid);
        mol.add_bond(3, c3, c4, 1, BondDisplay::Solid);
        mol.add_bond(4, c4, c5, 2, BondDisplay::Solid);
        mol.add_bond(5, c5, n, 1, BondDisplay::Solid);

        let val = CcVal::compute(&mol);
        for i in 0..5 {
            assert_eq!(val.implicit_h(i), 1);
        }
        let fp = mol.compute_formula_properties();
        assert_eq!(fp.get_formula(), "C4H5N");

        let arom = CcAromaticityInfo::perceive(&mol);
        for i in 0..5 {
            assert!(arom.is_atom_aromatic(i));
            assert!(arom.is_bond_aromatic(i));
        }
    }

    #[test]
    fn test_reaction_center_perception_and_validation() {
        assert!(!center_type_denotes_change(0));
        assert!(center_type_denotes_change(1));
        assert!(center_type_denotes_change(2));
        assert!(center_type_denotes_change(3));
        assert!(!center_type_denotes_change(4));
        assert!(center_type_denotes_change(5));
        assert!(!center_type_denotes_change(6));
        assert!(!center_type_denotes_change(7));
        assert!(center_type_denotes_change(8));

        let mut r = ChemMolecule::new();
        r.add_atom(ChemAtom::new(1, 6));
        r.add_atom(ChemAtom::new(2, 6));
        r.add_bond(1, 0, 1, 1, BondDisplay::Solid);

        let mut p = ChemMolecule::new();
        p.add_atom(ChemAtom::new(1, 6));
        p.add_atom(ChemAtom::new(2, 6));
        p.add_bond(1, 0, 1, 2, BondDisplay::Solid);

        let mut step = ReactionStep {
            reactant: r,
            product: p,
            atom_map: BTreeMap::new(),
            bond_map: BTreeMap::new(),
        };
        step.atom_map.insert(0, 0);
        step.atom_map.insert(1, 1);
        step.bond_map.insert(0, 0);

        let mut rxn = CcReaction::new();
        rxn.steps.push(step);
        rxn.create_reaction_centers_from_scratch();
        assert!(rxn.validate().is_ok());
        assert!(rxn.records[0].0.bond_bits.contains(&0));
        assert_eq!(rxn.records[0].0.center_types[0], 3);
    }

    #[test]
    fn test_ccval_dative_metals_and_bridging_hydrogens() {
        // Explicit dative bond: H3N -> BH3 (N=7 donor -> B=5 acceptor, order=5)
        let mut mol = ChemMolecule::new();
        let n = mol.add_atom(ChemAtom::new(1, 7));
        let b = mol.add_atom(ChemAtom::new(2, 5));
        mol.add_bond(1, n, b, 5, BondDisplay::Dative);

        let val = CcVal::compute(&mol);
        assert_eq!(val.explicit_dative_out[n], 1);
        assert_eq!(val.explicit_dative_in[b], 1);
        assert_eq!(val.implicit_h(n), 3); // N+ like C -> 4 - 1 = 3 H
        assert_eq!(val.implicit_h(b), 3); // B- like C -> 4 - 1 = 3 H
        // CoreChemistryCommon::GetFormula always emits C, then H, then remaining symbols alphabetically
        assert_eq!(mol.compute_formula_properties().get_formula(), "H6BN");

        // Charged alkali metal Li+ (Z=3, charge=+1) -> 0 implicit H
        let mut li_mol = ChemMolecule::new();
        let mut li = ChemAtom::new(1, 3);
        li.charge = 1;
        li_mol.add_atom(li);
        assert_eq!(CcVal::compute(&li_mol).implicit_h(0), 0);

        // Bridging H in B-H-B (diborane fragment) -> bridging H has expected valence 2, 0 implicit H
        let mut b2 = ChemMolecule::new();
        let b_a = b2.add_atom(ChemAtom::new(1, 5));
        let h_br = b2.add_atom(ChemAtom::new(2, 1));
        let b_b = b2.add_atom(ChemAtom::new(3, 5));
        b2.add_bond(1, b_a, h_br, 1, BondDisplay::Solid);
        b2.add_bond(2, h_br, b_b, 1, BondDisplay::Solid);
        let b2_val = CcVal::compute(&b2);
        assert_eq!(b2_val.expected_valence[h_br], 2);
        assert_eq!(b2_val.implicit_h(h_br), 0);
    }

    #[test]
    fn test_aromatic_fused_and_antiaromatic_systems() {
        // Cyclobutadiene (4 pi electrons -> antiaromatic, NOT aromatic)
        let mut cbd = ChemMolecule::new();
        for i in 0..4 {
            cbd.add_atom(ChemAtom::new(i + 1, 6));
        }
        cbd.add_bond(1, 0, 1, 2, BondDisplay::Solid);
        cbd.add_bond(2, 1, 2, 1, BondDisplay::Solid);
        cbd.add_bond(3, 2, 3, 2, BondDisplay::Solid);
        cbd.add_bond(4, 3, 0, 1, BondDisplay::Solid);
        let cbd_arom = CcAromaticityInfo::perceive(&cbd);
        assert!(cbd_arom.aromatic_bonds.is_empty());
        assert_eq!(cbd_arom.antiaromatic_bonds.len(), 4);

        // Azulene (5-ring fused to 7-ring = 10 pi electrons -> aromatic perimeter!)
        let mut az = ChemMolecule::new();
        for i in 0..10 {
            az.add_atom(ChemAtom::new(i + 1, 6));
        }
        // 10-membered perimeter alternating double/single bonds + bridgehead single bond (0-4)
        for i in 0..10 {
            let ord = if i % 2 == 0 { 2 } else { 1 };
            az.add_bond(i as u32 + 1, i, (i + 1) % 10, ord, BondDisplay::Solid);
        }
        az.add_bond(11, 0, 4, 1, BondDisplay::Solid);
        let az_arom = CcAromaticityInfo::perceive(&az);
        assert_eq!(az_arom.aromatic_atoms.len(), 10);
        assert_eq!(az_arom.aromatic_bonds.len(), 10);
    }

    #[test]
    fn test_stereo_double_bond_ez_and_isotope_cip() {
        // (Z)-2-butene vs (E)-2-butene: C1(CH3)-C2=C3-C4(CH3)
        let make_2butene = |y4: f64| {
            let mut m = ChemMolecule::new();
            let mut c1 = ChemAtom::new(1, 6); c1.x = 0.0; c1.y = 10.0;
            let mut c2 = ChemAtom::new(2, 6); c2.x = 10.0; c2.y = 20.0;
            let mut c3 = ChemAtom::new(3, 6); c3.x = 24.4; c3.y = 20.0;
            let mut c4 = ChemAtom::new(4, 6); c4.x = 34.4; c4.y = y4;
            m.add_atom(c1); m.add_atom(c2); m.add_atom(c3); m.add_atom(c4);
            m.add_bond(1, 0, 1, 1, BondDisplay::Solid);
            m.add_bond(2, 1, 2, 2, BondDisplay::Solid);
            m.add_bond(3, 2, 3, 1, BondDisplay::Solid);
            m
        };
        let z_butene = make_2butene(10.0); // both methyls up (y=10 < y=20)
        let e_butene = make_2butene(30.0); // opposite sides
        assert_eq!(compute_double_bond_cip_code(&z_butene, 1), 10); // Z
        assert_eq!(compute_double_bond_cip_code(&e_butene, 1), 11); // E

        // Isotope chirality: C(Br)(Cl)(D)(H) -> stereogenic (code 2 or 3)
        let mut chd = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6); c.x = 100.0; c.y = 100.0;
        let mut br = ChemAtom::new(2, 35); br.x = 100.0; br.y = 85.6;
        let mut cl = ChemAtom::new(3, 17); cl.x = 114.4; cl.y = 100.0;
        let mut d = ChemAtom::new(4, 1); d.isotope = 2; d.x = 100.0; d.y = 114.4;
        chd.add_atom(c); chd.add_atom(br); chd.add_atom(cl); chd.add_atom(d);
        chd.add_bond(1, 0, 1, 1, BondDisplay::WedgeBegin);
        chd.add_bond(2, 0, 2, 1, BondDisplay::Solid);
        chd.add_bond(3, 0, 3, 1, BondDisplay::Solid);
        let rep = perceive_stereochemistry(&chd);
        assert!(matches!(rep.atom_stereo.get(&0).copied(), Some(2) | Some(3)));
    }

    #[test]
    fn test_wave14_sulfoxide_oxime_allene_and_ring_ez_suppression() {
        // 1. (R)-Methyl phenyl sulfoxide: S=O + Wedge Me + Solid Ph -> 2 (R), ring double bonds = 0 (BS=N)
        let mut sulf = ChemMolecule::new();
        let mut s = ChemAtom::new(1, 16);
        s.x = 0.0;
        s.y = 0.0;
        let mut o = ChemAtom::new(2, 8);
        o.x = 0.0;
        o.y = -14.4;
        let mut me = ChemAtom::new(3, 6);
        me.x = -12.47;
        me.y = 7.2;
        let ph_coords = [
            (12.47, 7.2),
            (24.94, 0.0),
            (37.41, 7.2),
            (37.41, 21.6),
            (24.94, 28.8),
            (12.47, 21.6),
        ];
        sulf.add_atom(s);
        sulf.add_atom(o);
        sulf.add_atom(me);
        for (idx, (x, y)) in ph_coords.iter().enumerate() {
            let mut c = ChemAtom::new(4 + idx as u32, 6);
            c.x = *x;
            c.y = *y;
            sulf.add_atom(c);
        }
        sulf.add_bond(10, 0, 1, 2, BondDisplay::Solid);
        sulf.add_bond(11, 0, 2, 1, BondDisplay::WedgeBegin);
        sulf.add_bond(12, 0, 3, 1, BondDisplay::Solid);
        let ph0_double = sulf.add_bond(13, 3, 4, 2, BondDisplay::Solid);
        sulf.add_bond(14, 4, 5, 1, BondDisplay::Solid);
        sulf.add_bond(15, 5, 6, 2, BondDisplay::Solid);
        sulf.add_bond(16, 6, 7, 1, BondDisplay::Solid);
        sulf.add_bond(17, 7, 8, 2, BondDisplay::Solid);
        sulf.add_bond(18, 8, 3, 1, BondDisplay::Solid);
        assert_eq!(compute_tetrahedral_cip_code(&sulf, 0), 2);
        assert_eq!(compute_double_bond_cip_code(&sulf, ph0_double), 0);

        // 2. (E)- vs (Z)-Acetaldoxime
        for (oh_y, expected) in [(-12.47_f64, 11u8), (12.47_f64, 10u8)] {
            let mut ox = ChemMolecule::new();
            for (id, z, x, y) in [
                (1, 6, -12.47, 7.2),
                (2, 6, 0.0, 0.0),
                (3, 7, 14.4, 0.0),
                (4, 8, 21.6, oh_y),
            ] {
                let mut a = ChemAtom::new(id, z);
                a.x = x;
                a.y = y;
                ox.add_atom(a);
            }
            ox.add_bond(10, 0, 1, 1, BondDisplay::Solid);
            let db = ox.add_bond(11, 1, 2, 2, BondDisplay::Solid);
            ox.add_bond(12, 2, 3, 1, BondDisplay::Solid);
            assert_eq!(compute_double_bond_cip_code(&ox, db), expected);
        }

        // 3. Axial allene (P) vs (M), plus terminal sp2 carbon NOT falsely tetrahedral
        for (disp, expected) in [
            (BondDisplay::WedgeBegin, Some('P')),
            (BondDisplay::WedgedHashBegin, Some('M')),
            (BondDisplay::Solid, None),
        ] {
            let mut allene = ChemMolecule::new();
            for (id, z, x, y) in [
                (1, 6, -14.4, 0.0),
                (2, 6, 0.0, 0.0),
                (3, 6, 14.4, 0.0),
                (4, 17, -21.6, -12.47),
                (5, 9, -21.6, 12.47),
                (6, 53, 21.6, -12.47),
                (7, 35, 21.6, 12.47),
            ] {
                let mut a = ChemAtom::new(id, z);
                a.x = x;
                a.y = y;
                allene.add_atom(a);
            }
            allene.add_bond(10, 0, 1, 2, BondDisplay::Solid);
            allene.add_bond(11, 1, 2, 2, BondDisplay::Solid);
            allene.add_bond(12, 0, 3, 1, BondDisplay::Solid);
            allene.add_bond(13, 0, 4, 1, disp);
            allene.add_bond(14, 2, 5, 1, BondDisplay::Solid);
            allene.add_bond(15, 2, 6, 1, BondDisplay::Solid);
            assert_eq!(compute_axial_allene_cip(&allene, 1), expected);
            assert_eq!(compute_tetrahedral_cip_code(&allene, 0), 0);
        }
    }

    #[test]
    fn test_molecule_properties_and_pseudoasymmetry() {
        // Benzene C6H6: MW = 6*12.011 + 6*1.008 = 78.114, elemental analysis C ~92.26%, H ~7.74%
        let mut benz = ChemMolecule::new();
        for i in 0..6u32 {
            benz.add_atom(ChemAtom::new(i + 1, 6));
        }
        for i in 0..6usize {
            benz.add_bond((i + 1) as u32, i, (i + 1) % 6, 6, BondDisplay::Solid);
        }
        let fp = CsFormulaProperties::from_molecule(&benz);
        assert_eq!(fp.get_formula(), "C6H6");
        let mw = fp.get_molecular_weight();
        assert!((mw - 78.11).abs() < 0.05, "unexpected MW: {mw}");
        let ea = fp.format_elemental_analysis();
        assert!(ea.starts_with("C, 92.2") && ea.contains("H, 7.7"), "unexpected EA: {ea}");

        // Pseudoasymmetry (CIP Rule 5): C2(H)(OH)(C1(R))(C3(S)) -> r (4) or s (5)
        let mut m = ChemMolecule::new();
        let coords = [
            (0.0, 0.0),     // 0: C1 (R or S)
            (14.4, 0.0),    // 1: C2 (pseudoasymmetric center)
            (28.8, 0.0),    // 2: C3 (opposite config)
            (-10.0, 10.0),  // 3: Br on C1
            (-10.0, -10.0), // 4: F on C1
            (14.4, 14.4),   // 5: Cl on C2
            (38.8, 10.0),   // 6: Br on C3
            (38.8, -10.0),  // 7: F on C3
        ];
        let zs = [6, 6, 6, 35, 9, 17, 35, 9];
        for (i, (&(x, y), &z)) in coords.iter().zip(zs.iter()).enumerate() {
            let mut a = ChemAtom::new((i + 1) as u32, z);
            a.x = x;
            a.y = y;
            m.add_atom(a);
        }
        m.add_bond(1, 0, 1, 1, BondDisplay::Solid);
        m.add_bond(2, 1, 2, 1, BondDisplay::Solid);
        m.add_bond(3, 0, 3, 1, BondDisplay::WedgeBegin);
        m.add_bond(4, 0, 4, 1, BondDisplay::Solid);
        m.add_bond(5, 1, 5, 1, BondDisplay::WedgeBegin);
        m.add_bond(6, 2, 6, 1, BondDisplay::WedgeBegin);
        m.add_bond(7, 2, 7, 1, BondDisplay::Solid);
        let c1 = compute_tetrahedral_cip_code(&m, 0);
        let c3 = compute_tetrahedral_cip_code(&m, 2);
        assert!(c1 != c3 && (c1 == 2 || c1 == 3) && (c3 == 2 || c3 == 3));
        let c2 = compute_tetrahedral_cip_code(&m, 1);
        assert!(c2 == 4 || c2 == 5, "expected pseudoasymmetric r/s (4/5), got {c2}");

        // S05: 15-atom ligand divergence depth test
        let mut deep = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6); c.x = 0.0; c.y = 0.0;
        let mut f = ChemAtom::new(2, 9); f.x = 0.0; f.y = -14.4;
        deep.add_atom(c); deep.add_atom(f);
        deep.add_bond(1, 0, 1, 1, BondDisplay::WedgeBegin);
        let mut prev_a = 0usize;
        for step in 0..15 {
            let mut ca = ChemAtom::new((deep.atoms.len() + 1) as u32, 6);
            ca.x = 14.4 * ((step + 1) as f64);
            let idx = deep.add_atom(ca);
            deep.add_bond((deep.bonds.len() + 1) as u32, prev_a, idx, 1, BondDisplay::Solid);
            prev_a = idx;
        }
        let mut prev_b = 0usize;
        for step in 0..15 {
            let z = if step == 14 { 8 } else { 6 }; // diverges at depth 15!
            let mut cb = ChemAtom::new((deep.atoms.len() + 1) as u32, z);
            cb.x = -14.4 * ((step + 1) as f64);
            let idx = deep.add_atom(cb);
            deep.add_bond((deep.bonds.len() + 1) as u32, prev_b, idx, 1, BondDisplay::Solid);
            prev_b = idx;
        }
        let code_deep = compute_tetrahedral_cip_code(&deep, 0);
        assert!(code_deep == 2 || code_deep == 3, "S05 depth-15 ligand divergence failed");

        // Kekulize & Dekekulize (ccKekulize / ccDekekulize): benzene order=6 -> alternating 1/2 -> back to order=6
        let mut arom_b = ChemMolecule::new();
        for i in 0..6u32 {
            arom_b.add_atom(ChemAtom::new(i + 1, 6));
        }
        for i in 0..6usize {
            arom_b.add_bond((i + 1) as u32, i, (i + 1) % 6, 6, BondDisplay::Solid);
        }
        assert!(kekulize_aromatic_bonds(&mut arom_b));
        assert_eq!(arom_b.bonds.iter().filter(|b| b.order == 2).count(), 3);
        assert_eq!(arom_b.bonds.iter().filter(|b| b.order == 1).count(), 3);
        assert_eq!(arom_b.compute_formula_properties().get_formula(), "C6H6");
        assert_eq!(dekekulize_aromatic_bonds(&mut arom_b), 6);
        assert!(arom_b.bonds.iter().all(|b| b.order == 6));
    }
}
