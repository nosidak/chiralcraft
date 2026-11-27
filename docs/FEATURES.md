# Complete Feature Reference

ChiralCraft unifies publication-grade 2D chemical drawing, real-time cheminformatics, medicinal chemistry, multi-spectroscopy, retrosynthesis, organometallics, biopolymers, 3D mechanics, crystallography, and quantum chemistry into one native pure-Rust desktop application.

---

## 1. Publication-Grade 2D Drawing & Visual Polish

### IUPAC 2008 Rendering Engine
- **Smart Terminal Heteroatom Inversion (`GR-2.1`)**: Automatically flips left-facing terminal groups (`HO—`, `H₂N—`, `HS—`, `H₃C—`) so bonds always attach to the chemically bonded atom.
- **Typographic Formulas**: Proper superscript formal charges (`+`, `−`, `2+`, `3+`), superscript atomic mass isotopes (`¹³C`, `²H`), and subscript implicit hydrogen counts.
- **Linear & Ring Geometry**: Automatic `180°` extension for sp-hybridized alkynes/nitriles, `120°` chain angles for sp²/sp³ chains, ring-interior shortened double bonds, and toggleable aromatic circles (`Alt+K`).
- **Crossing-Bond Bridge Gaps**: Automatically detects overlapping bond pairs (such as in adamantane, norbornane, cubane, and macrocycles) and inserts configurable white/background clearance gaps behind front bonds.
- **3D Rear-Bond Depth Cueing**: Configurable rear-atom/bond opacity (`0.15–1.0`) for 3D projections and bridged cages.

### 16 Built-In Journal Stylesheets
Apply any stylesheet to a new document (`File › Open Style Sheets`), existing document (`File › Apply Document Settings from`), or active structure (`Object › Apply Object Settings from`):
- `ACS Document 1996` (`14.40 pt`)
- `Nature Research Document` (`12.00 pt`)
- `RSC (1 column)` & `RSC (2 column)`
- `Wiley Document`
- `JOC (Journal of Organic Chemistry)`
- `Angewandte Chemie`
- `Science of Synthesis`
- `Phytomedicine Document`
- `J. Het. Chem.`
- `J. Mol. Mod. (1 column)` & `J. Mol. Mod. (2 column)`
- `Dark Mode`
- `New Slide` (presentation bold)
- `New A4 Document`
- `Standard 30 pt Document`

### 11 Built-In Clean-Room Template Libraries & Special Cores
- **Chemical Libraries**: Aromatics & Heterocycles · Amino Acids · Carbohydrates · Nucleobases · Pharmaceuticals · Functional Groups · Steroids & Polycycles · Haworth & Macrocycles · Organometallic Catalysts.
- **Vector Graphics Libraries**: Lab Glassware & Apparatus (flasks, condensers, Schlenk lines, columns) · BioGraphics & Lab Shapes (membranes, antibodies, helices, vessels).
- **Special Core Builders (`Structure › Insert Special Core`)**: Haworth β-D-Glucopyranose, Haworth β-D-Ribofuranose, Porphine Macrocycle, 18-Crown-6, 3D-projected Adamantane, and Norbornane.
- **Catalytic Cycles (`Structure › Catalytic Cycle`)**: Arrange any 2+ selected molecules into a circular or elliptical catalytic cycle with curved connecting arrows, or stamp a complete 4-step Pd starter catalytic cycle.

### Publication TOC Composer, TLC Plates & Data Tables
- **TOC Boundary Frames**: Exact journal Graphical Abstract / Table-of-Contents bounding frames (`ACS 3.25 × 1.75 in`, `Nature`, `RSC`, `Angewandte`) with live overflow detection.
- **Smart Arrow Conditions & Serpentine Wrapping**: Stack reagent/condition captions centered above and below any reaction arrow; wrap multi-step reaction sequences into tidy serpentine rows.
- **Interactive TLC Plate Editor**: Configure eluent label, origin baseline fraction, spot `Rf` positions, lane labels, and custom RGB spot colors — persisted losslessly in CDXML and SVG.
- **Editable On-Canvas Tables**: Multi-row/column data tables (e.g. reaction optimization tables) editable directly in the Object Inspector and exported cleanly to SVG.

---

## 2. Real-Time Cheminformatics, Naming & Databases

- **Chemistry Analysis & Properties**: Formula (Hill system), Exact Mass, Molecular Weight, Elemental Analysis (`C, H, N, O...`), Ertl tPSA, Wildman–Crippen LogP & Molar Refractivity, LogS solubility, boiling/melting point estimates, and pKa estimates.
- **Live Linked Captions**: Attach `Formula`, `Molecular Weight`, `Exact Mass`, `IUPAC Name`, or `Full Summary` captions that dynamically recompute whenever you edit atoms or bonds.
- **Chemical Structure Warnings (`Structure › Check Structure`)**: 13 configurable validation checks covering valence violations, charge mismatches, isolated fragments, overlapping atoms, ambiguous stereocenters, and linear geometry errors.
- **Stereochemistry**: Automatic CIP `(R)`/`(S)` and `(E)`/`(Z)` perception, ring `(E)`/`(Z)` suppression, sulfoxide/oxime/allene/atropisomer handling, and MDL V3000 **Enhanced Stereochemistry** (`abs`, `or1`, `&1`).
- **Hybrid Structure ↔ Name**: Fast pure-Rust offline IUPAC naming + automatic live **PubChem** IUPAC & common drug name lookup (`Alt+Ctrl+N`) and Name-to-Structure (`Shift+Ctrl+N`), plus user-defined custom dictionary entries.
- **Live PubChem & ChEMBL Search (`Shift+Ctrl+F`)**: Search by name, CAS, CID, ChEMBL ID, formula, exact SMILES, 2D Tanimoto similarity, or substructure; inspect clinical phase and drug metadata; import clean 2D MOL/SDF structures normalized to your active stylesheet.
- **Chemistry-Aware Document Diff (`File › Compare Documents...`)**: Compare any two open tabs or files chemically — detecting added, removed, and modified structures, formula/mass deltas, and stereochemical changes.

---

## 3. Medicinal Chemistry, ADMET, PAINS & HRMS

- **MedChem & ADMET Scorecard**: Live calculation of MW, LogP, tPSA, HBD, HBA, Rotatable Bonds, Heavy Atoms, **Fsp³**, **Lipinski Rule of 5**, **Veber**, **Ghose**, **Pfizer CNS MPO (0–6)**, **QED (0–1)**, and **Synthetic Accessibility (SA, 1–10)**.
- **PAINS & Brenk Structural Alerts**: Flags Michael acceptors, epoxides, nitroaromatics, acyl halides, rhodanines, quinones, hydrazines, aldehydes, peroxides, and other reactive/pan-assay interference motifs with live red/amber atom highlights on canvas.
- **HRMS Isotope Envelope & Adduct Simulator**: Exact isotopic convolution across all 118 IUPAC elements (`³⁵Cl/³⁷Cl`, `⁷⁹Br/⁸¹Br`, `¹⁰B/¹¹B`, `³²S/³⁴S`, `Pd`, `Pt`, etc.), interactive resolution control (`R = 5,000–100,000`), Gaussian profile + centroid stick display, 10 ESI+/ESI− adducts, and observed `m/z` ppm error calculator.
- **Virtual Screening Library Browser**: Load `.sdf`, `.smi`, or `.csv` libraries (or built-in demo library), filter/sort by Lipinski/MPO/QED/PAINS, and stamp hits onto the canvas.

---

## 4. Multi-Spectroscopy Studio

Open via **`Structure › Multi-Spectroscopy Studio (FT-IR / UV–Vis / 2D NMR)...`** or **Right Studio Dock → `Props` / `Spectra` → `FT-IR / Spectra`**.

- **¹H-NMR (400 MHz) & ¹³C-NMR (100 MHz)**: Interactive Lorentzian spectrum viewer with reversed chemical-shift axis, multiplicity & coupling estimates, clickable peak ↔ atom canvas highlighting, and collision-free radial on-canvas `δ (ppm)` shift badges (`Structure › NMR Shift Labels`).
- **Experimental JCAMP-DX (`.jdx`) & CSV Overlay**: Browse and load real spectrometer `.jdx`/`.dx`/`.csv`/`.dpt` files (`XYDATA`, IUPAC ASDF/SQZ/DIF/DUP, `PEAK TABLE`) — or fetch directly from **NIST Chemistry WebBook** — with live **Hit Quality Index (HQI %)** and Pearson `r` match scoring.
- **Multi-Backend FT-IR Studio (`4000–400 cm⁻¹`)**:
  - **Instant Rules Engine**: Ring-strain, conjugation, H-bonding, and substitution-aware vibrational band prediction with per-band FWHM broadening.
  - **Local Quantum Harmonic IR (`GFN2-xTB` / `GFN-FF` / `ORCA` / `Gaussian`)**: Run local `xtb --ohess` or import `.out`, `.log`, `g98.out`, `vibspectrum`, or `.molden` outputs with anharmonic scaling presets (`0.967 B3LYP`, `0.970 GFN2-xTB`, `0.975 r²SCAN-3c`, `0.950 ωB97X-D`), Lorentzian/Gaussian/Pseudo-Voigt lineshapes, scroll-wheel zoom, and **3D normal-mode vibration animation** (`Play 3D`).
  - **Publication Export**: Save vector `.svg`, `.jdx`, or `.csv` spectra via native OS file dialogs, or paste formatted ACS IR strings onto the canvas.
- **UV–Vis Chromophore Predictor**: Conjugation-aware `λmax` (nm), molar absorptivity `ε`, electronic transition type (`π→π*`, `n→π*`), and perceived solution color swatch.
- **2D NMR Contour Maps (`COSY`, `HSQC`, `HMBC`)**: Interactive 2D correlation plots where clicking any cross-peak highlights both coupled atoms on the canvas.

---

## 5. Synthesis Planning, Retrosynthesis, Isomers, Kinetics & ELN Safety

### How to Use the Retrosynthetic Disconnection Advisor

1. **Open the Advisor** using any of three quick paths:
   - Top menu: **`Structure › Retrosynthetic Disconnection Advisor...`**
   - Right Studio Dock: **`Props` tab → `Retrosynthesis`** (always visible, even on a blank canvas)
   - Command Palette: press **`:`** or **`F4`**, type **`retro`**, and press `Enter`
2. **Choose or Draw a Target Molecule**:
   - Draw any molecule on the canvas (or paste a SMILES/name via `Shift+Ctrl+N` / `Alt+Ctrl+P`) — disconnections run **automatically** when you open the advisor or edit the molecule.
   - Want to test it right away on a blank canvas? Click one of the built-in **1-click examples** inside the advisor window: **`Acetanilide`**, **`Aspirin`**, **`Biphenyl`**, or **`Sulfanilamide`**.
3. **Inspect & Stamp a Synthetic Route**:
   - Click any ranked disconnection row to highlight the **strategic target bond in dashed red** directly on the canvas and inspect the suggested **forward reagents** and **synthons**.
   - Click **`Stamp Scheme (⇒)`** to place a publication-ready retrosynthesis scheme (`Target ⇒ Precursor 1 + Precursor 2` with forward reagents) onto the canvas.
4. **Supported Strategic Disconnection Classes (9 Built-In Transforms)**:
   - **Amide Coupling** (HATU / EDC·HCl, DIPEA)
   - **Esterification** (EDC/DMAP or Fischer)
   - **Suzuki–Miyaura Biaryl Cross-Coupling** (`Pd(PPh₃)₄`, `K₂CO₃`)
   - **Buchwald–Hartwig C–N Amination** (`Pd₂(dba)₃`, XPhos, `NaOtBu`)
   - **Sulfonamide Coupling** (`Et₃N`, DCM)
   - **Williamson Ether Synthesis** (`K₂CO₃` or `NaH`)
   - **Reductive Amination** (`NaBH(OAc)₃`, AcOH)
   - **Friedel–Crafts Acylation** (`AlCl₃`, DCM)
   - **Wittig Olefination** (`Ph₃P=CHR`, THF)
5. **CLI & MCP Access**: Run `chiralcraft retro target.mol` (or `chiralcraft retro "CC(=O)Nc1ccccc1"`) from the terminal, or call `retrosynthesis_suggest` via `chiralcraft --mcp`.

### Additional Synthesis & ELN Tools
- **Markush R-Group Enumerator & SAR Grid**: Enumerate `R1`/`R2` substituent combinations onto any scaffold, sort the generated SAR library by MW/LogP/MPO/QED/Alerts, stamp visual grids onto the canvas, or export to `.sdf`/`.csv`.
- **pH Microspecies, Tautomers, Stereoisomers & Resonance**:
  - Live `pH 1.0–14.0` protonation state explorer with calculated isoelectric point (`pI`).
  - Prototropic tautomer enumerator (keto–enol, imine–enamine, amide–imidic acid, lactam–lactim).
  - Full `2ⁿ` stereoisomer enumerator with `(R)`/`(S)` labels.
  - Ranked resonance contributors & interactive curved-arrow electron-pushing mechanistic stepper.
- **Reaction Energy Profiles (`ΔG‡`) & Eyring/Boltzmann Kinetics**: Interactive free-energy diagram editor with optional competing Pathway B curves, Eyring transition-state rate constant `k` & half-life `t₁/₂` at arbitrary temperature, Boltzmann selectivity (`% ee` / `er` from `ΔΔG‡`), and lossless SVG/CDXML stamping.
- **Stoichiometry Grid, AutoNumbering & ELN Supporting Info**:
  - Automatic reactant/product detection, mmol/mg/equiv recalculation, and CSV export.
  - Searchable **30-entry Lab Solvent & Reagent GHS Safety Handbook** (MW, density, boiling point, hazards, automatic liquid mL volume calculator).
  - One-click ACS/JOC Supporting Information experimental procedure + formatted ¹H-NMR, ¹³C-NMR, FT-IR, and HRMS characterization paragraph.

---

## 6. Organometallic, Biopolymer & Polymer Chemistry

- **Organometallic 18-Electron & Green CBC (`MLxXy`) Engine**: Computes formal oxidation state, `dⁿ` count, L-type/X-type/η⁵/κ² ligand breakdown, metal–metal bonds, total valence electrons, saturation status (`18e⁻`, `16e⁻ d⁸ square planar`, coordinatively unsaturated), and predicted geometry.
- **Biopolymer Studio (`HELM` / `FASTA`)**: Sequence builder for Peptides, DNA, and RNA, monomer atom stamping, serpentine biopolymer layout cleanup, and editable **Gel Electrophoresis Plates (GEP)**.
- **Polymer `Mn`/`Mw`/`Đ` & BigSMILES**: Repeat-unit molecular weight, dispersity (`Đ = Mw/Mn`), degree of polymerization, and BigSMILES string generator.

---

## 7. Interactive 3D, Crystallography & Quantum Chemistry

- **Interactive 3D Molecular & Crystal Viewer**: Ball & Stick, Stick, Space-Filling (CPK), and Wireframe modes, drag-to-orbit, zoom slider, dark/light stage toggle, crystallographic `a`/`b`/`c`/`Iso` axis snap, and **Snapshot to 2D Canvas**.
- **UFF-Lite 3D Force-Field Optimizer & Measurement**: Pure-Rust bond-stretch, angle-bend, torsion, and van der Waals minimization plus click-to-measure **bond distances (Å)**, **bond angles (°)**, and **dihedral angles (°)**.
- **Quantum Chemistry Bridge (ORCA / Gaussian / xTB)**: Generate ready-to-run ORCA `.inp`, Gaussian `.gjf`, and `GFN2-xTB` input decks (Opt+Freq, Single Point, TD-DFT, GIAO NMR, implicit solvents, charge/multiplicity), and parse `.out`/`.log` output logs for SCF energy, HOMO–LUMO gap, dipole moment, DFT IR frequencies, and optimized 3D coordinates.
