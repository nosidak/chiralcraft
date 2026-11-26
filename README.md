<p align="center">
  <img src="assets/chiralcraft-icon.svg" alt="ChiralCraft Icon" width="128" height="128" />
</p>

# ChiralCraft

**ChiralCraft** is a fast, native, pure-Rust desktop chemical structure editor and cheminformatics CLI built with `eframe`/`egui`.

## Highlights
- **Publication-Grade 2D Drawing**: Follows IUPAC 2008 (`GR-1`–`GR-10`) and Nature Research guidelines (automatic left-facing terminal heteroatom inversion `HO—`, `H₂N—`, `HS—`, superscript charges & isotopes, `180°` linear alkynes/nitriles, ring-interior double bonds, wedged/hashed/wavy stereobonds, automatic **crossing-bond bridge gaps**, and 16 built-in journal stylesheets including ACS 1996, Nature Research, RSC, Wiley, JOC, and Angewandte Chemie).
- **Catalytic Cycles & Special Cores**: One-click circular/elliptical catalytic cycle layout (`Structure > Catalytic Cycle`), starter 4-step catalytic cycle generator, plus **Haworth sugar projections** (pyranose/furanose), **porphine**, **18-crown-6**, **adamantane**, **norbornane**, and modular **lab glassware & gas-line apparatus** across **10 built-in clean-room template libraries**.
- **Real-Time Chemistry & Live Linked Captions**: Periodic & isotope tables, exact mass, molecular weight, elemental analysis, live valence/charge/stereo/geometry chemical warnings, aromaticity perception, CIP `(R)/(S)/(E)/(Z)` stereochemistry, enhanced stereochemistry (`abs`, `or1`, `&1`), and **Live Linked Property & Name Captions** that automatically stay in sync when you edit a molecule.
- **Interactive 3D & Crystal Viewer**: Drag-to-orbit Ball & Stick, Stick, Space-Filling (CPK), and Wireframe modes with VESTA-inspired coloring, quick 2D→3D conformer lift, native `.xyz`, `.cif`, `.pdb`, and VASP `POSCAR`/`CONTCAR` import with **crystallographic unit-cell wireframe** & `a`/`b`/`c` axis snap, **"Snapshot to 2D Canvas"**, and adjustable **3D Rear Opacity** depth cueing.
- **Command Palette (`:` / `F4`) & Smart Reverse Group Typing**: Instant fuzzy action search + parameterized commands (`rotate 30`, `scale 150`, `smiles CCO`, `name aspirin`, `style nature`, `theme slide`, `nmr stamp`, `3d`, `diff`, `dpi 600`), plus automatic resolution of left-facing reversed groups (`HO`, `H2N`, `O2N`, `MeO`, `AcO`, `F3C`, `BocHN`, `HO2C`).
- **Chemistry-Aware Document Diff & On-Canvas NMR Labels**: Compare tabs or files chemically (`File > Compare Documents...` or `chiralcraft diff a.svg b.svg`) highlighting added/removed/modified structures and stereo changes, and stamp predicted ¹H / ¹³C NMR `δ (ppm)` shift badges right onto the canvas (`chiralcraft nmr mol.mol`).
- **Editable SVG, PNG & Native Vector PDF Round-Trip**: Embeds the lossless CDXML chemical graph inside **SVG** (`<chiralcraft:document>`), **PNG** (`tEXt` chunk `ChiralCraftDocument`), and native pure-Rust **PDF 1.4** (`%CHIRALCRAFT_CDXML_B64:`), with configurable DPI (`150/300/600`), padding, and copy background preferences.
- **Interactive TLC Plate & Table Inspector**: Edit TLC eluent labels, baseline, spot $R_f$ values, lane labels, spot colors, and multi-cell tables live on canvas and in exported SVGs.
- **Independent Dark / Light UI & Dark-Slide Canvas Themes**: Switch freely between **Classic Light**, **Dark UI + White Canvas**, and **Dark UI + Dark Slide Canvas** (with automatic white-ink transparent SVG export for dark presentations).
- **Live PubChem & ChEMBL Search & Import**: Search compounds by name, CAS, formula, CID, ChEMBL ID, exact SMILES, 2D Tanimoto similarity, or substructure (`Shift+Ctrl+F`), preview 2D structures with clinical/drug metadata, and import curated 2D MOL/SDF structures scaled to your active journal stylesheet.
- **Hybrid Structure ↔ Name**: Instant pure-Rust offline IUPAC naming + automatic live PubChem IUPAC & common drug name resolution (`Alt+Ctrl+N`), with live name display in the status bar and Analysis panel.
- **Spectroscopy, Biopolymers & Reaction Tools**: Interactive ¹H-NMR & ¹³C-NMR prediction viewer with live peak↔atom highlighting, HELM/FASTA peptide/DNA/RNA builder, Reaction Stoichiometry Grid with CSV export, AutoNumbering, Reaction Atom Mapping, and Gel Electrophoresis Plate (GEP).
- **Medicinal Chemistry, ADMET & PAINS Scorecard**: Live Lipinski, Veber, Ghose, Pfizer CNS MPO (0–6), QED drug-likeness, Fsp³, Synthetic Accessibility (SA), and PAINS/Brenk structural alerts with real-time red/amber atom highlighting directly on the canvas.
- **HRMS Isotope Envelope & Adduct Simulator**: Exact natural-abundance isotopic convolution across all 118 IUPAC elements (`Cl`/`Br`/`S`/`B`/`Si`/`Pd`/`Pt` patterns), interactive Gaussian profile & centroid stick plot, ESI+/ESI- adduct table, and observed $m/z$ ppm-error calculator.
- **pH Microspecies, Tautomers, Stereoisomers & Resonance**: Physiological pH (`1.0–14.0`) protonation slider with isoelectric point ($\text{p}I$), prototropic tautomer enumerator, $2^n$ stereoisomer enumerator, ranked resonance contributors, and curved-arrow electron-pushing mechanistic stepper.
- **Markush R-Group Library Enumerator & Retrosynthesis Advisor**: Enumerate `R1`/`R2` scaffolds into sortable SAR grids (with SDF/CSV export), and analyze targets for strategic disconnections (Amide, Ester, Suzuki, Buchwald–Hartwig, Sulfonamide, Williamson, Reductive Amination, Friedel–Crafts, Wittig) with one-click `⇒` retrosynthesis scheme stamping.
- **Publication TOC Composer & JCAMP-DX Experimental NMR Overlay**: Exact journal TOC boundary frames (`ACS`, `Nature`, `RSC`, `Angewandte`) with live overflow checks, smart arrow condition stacking, serpentine multi-step scheme wrapping, and `.jdx` experimental NMR overlay behind predicted ¹H/¹³C traces.
- **Multi-Spectroscopy Studio (FT-IR, UV–Vis & 2D NMR)**: Simulated FT-IR (`4000–400 cm⁻¹` `%T`) with clickable bond-highlighting band assignments, UV–Vis $\lambda_{\max}$ & perceived solution color swatch, and interactive 2D NMR contour plots (**COSY**, **HSQC**, **HMBC**) highlighting both coupled atoms on canvas.
- **Organometallic 18-Electron Engine & Polymer / BigSMILES Studio**: Green CBC (`MLxXy`) transition-metal analyzer (oxidation state, $d^n$, 18e⁻/16e⁻ check, coordination geometry) + built-in **Organometallic Catalysts** library (#11: `Pd(PPh₃)₄`, Cisplatin, Wilkinson's, Vaska's, Ferrocene, Grubbs Ru), alongside a Polymer $M_n$/$M_w$/Đ & **BigSMILES** calculator.
- **Reaction Energy Profile (`ΔG‡`) & Eyring/Boltzmann Kinetics**: Interactive canvas & SVG free-energy diagram builder with competing Pathway B curves, Eyring rate constant $k$ & half-life $t_{1/2}$, and Boltzmann selectivity (`% ee` / `er` from `ΔΔG‡`).
- **3D UFF Molecular Mechanics Optimizer & Geometry Measurement**: Pure-Rust 3D force-field energy minimization inside the interactive 3D Viewer plus click-to-measure **bond distances (Å)**, **angles (°)**, and **dihedrals (°)**.
- **ELN Supporting Info Generator & GHS Lab Safety Handbook**: Searchable solvent/reagent safety & physical-property database with automatic liquid volume (mL) calculation and one-click ACS/JOC experimental procedure + NMR/IR/HRMS characterization text generation.
- **Native LaTeX (`\chemfig` / `\schemestart`) / TikZ Export & SDF Virtual Screening**: Export compilable LaTeX/TikZ code and screen multi-compound `.sdf`/`.smi`/`.csv` libraries by MW, LogP, CNS MPO, Lipinski, and PAINS alerts.
- **Computational Chemistry Bridge (ORCA / Gaussian / xTB)**: Export ready-to-run `.inp`, `.gjf`, and `GFN2-xTB` input decks straight from any 2D/3D structure and parse `.out`/`.log` results (SCF energy, HOMO–LUMO gap, dipole moment, DFT IR modes, and optimized 3D geometry).
- **Unified Right-Side Studio Dock & Desktop Integration**: Tabbed sidebar (`Props`, `MedChem`, `Spectra`, `DFT`, `Inspect`), persistent `~/.config/chiralcraft/config.json`, `File > Open Recent`, OS drag-and-drop, `Ctrl+Scroll` zoom-to-cursor, middle-mouse pan, and one-command Linux `.desktop` + MIME installer (`./scripts/install-linux-desktop.sh`).
- **22+ Chemical & Vector Formats + 17-Tool Built-In MCP Server**: Native read/write support for CDXML, CDX, MOL V2000/V3000, SDF, RXN V2000/V3000, SMILES, InChI, CML, CT, TGF, SKC, **XYZ**, **CIF**, **PDB**, **VASP POSCAR/CONTCAR**, **ORCA/Gaussian/xTB**, **JCAMP-DX (.jdx)**, **Editable SVG/PNG/PDF** (with full extension persistence for energy profiles, tables, and TLC plates), **LaTeX (.tex)**, EPS, EMF, and 3MF — plus a built-in stdio **Model Context Protocol (`--mcp`)** server with 17 tools for AI agents.

## Building, Installing & Running
```bash
cargo build -p chiralcraft-app --release
./scripts/install-linux-desktop.sh   # installs ~/.local/bin/chiralcraft, icon, MIME types & launcher
chiralcraft
```

### CLI & MCP Server Mode
```bash
./target/release/chiralcraft analyze molecule.mol
./target/release/chiralcraft medchem molecule.mol
./target/release/chiralcraft hrms molecule.mol "[M+H]+"
./target/release/chiralcraft ir molecule.mol
./target/release/chiralcraft 18e complex.cdxml
./target/release/chiralcraft retro molecule.mol
./target/release/chiralcraft nmr molecule.mol [--13c]
./target/release/chiralcraft quantum-in molecule.mol orca
./target/release/chiralcraft quantum-parse job.out
./target/release/chiralcraft export-tex scheme.cdxml
./target/release/chiralcraft screen library.sdf
./target/release/chiralcraft diff old.svg new.svg
./target/release/chiralcraft convert input.cif output.pdf
./target/release/chiralcraft inspect-media diagram.png
./target/release/chiralcraft --edit-svg /tmp/embedded-diagram.svg
./target/release/chiralcraft --mcp
```

See [docs/WORDCRAFT_DECKCRAFT_INTEGRATION.md](docs/WORDCRAFT_DECKCRAFT_INTEGRATION.md) for the WordCraft / DeckCraft round-trip protocol and [CHANGELOG.md](CHANGELOG.md) for release notes.

## Round-Trip Editable SVG Integration (WordCraft / DeckCraft / Office Suites)
Every SVG exported or copied (`Ctrl+C` or `Edit > Copy As > Editable SVG`) by ChiralCraft embeds the full lossless CDXML chemical graph inside standard SVG `<metadata>`:
```xml
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:chiralcraft="https://chiralcraft.dev/ns/1.0"
     data-chiralcraft-version="1.0" ...>
  <metadata>
    <chiralcraft:document format="cdxml" encoding="base64">...</chiralcraft:document>
  </metadata>
  <!-- publication-grade IUPAC vector paths & text -->
</svg>
```
- **Copy → Paste (`Ctrl+C` → `Ctrl+V`)**: `Edit > Copy` pushes `image/svg+xml` directly onto the Wayland (`wl-copy`) or X11 (`xclip`) clipboard with a transparent canvas background so structures paste cleanly into WordCraft documents and DeckCraft slides while remaining valid `.docx`/`.pptx` SVGs.
- **Double-Click Round-Trip Editing**: When a host app (WordCraft, DeckCraft, or a Linux office plugin) detects `<chiralcraft:document` inside a double-clicked SVG:
  1. Write the embedded SVG bytes to a temp path (e.g. `/tmp/wordcraft-chem-1.svg`).
  2. Launch `chiralcraft --edit-svg /tmp/wordcraft-chem-1.svg`.
  3. When the user saves (`Ctrl+S` / `Save & Return`) or closes ChiralCraft, ChiralCraft atomically writes the updated editable SVG back to that path — the host app reloads the file and updates the embedded graphic in place.

## License
GNU Affero General Public License v3.0 or later (`AGPL-3.0-or-later`) — see [LICENSE](LICENSE).

Copyright (c) 2026 [nosidak](https://github.com/nosidak).
