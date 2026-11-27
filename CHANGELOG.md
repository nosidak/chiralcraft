# Changelog

All notable changes to **ChiralCraft** are documented in this file.

## [1.1.0] - 2026-11-27

### Added
- **Multi-Backend FT-IR Studio**: Unified `IrPrediction` engine supporting Instant Rules (expanded with ring-strained carbonyls, anhydrides, aldehyde Fermi doublets, primary/secondary `NH₂` doublets, terminal alkynes, benzene OOP patterns, alkene OOP wags, `P=O`, `S–H`, `Si–O`, and `M–C≡O`), local async **GFN2-xTB / GFN-FF** (`--ohess`), **ORCA**, and **Gaussian** harmonic frequency runners with SMILES/method caching and safe GNOME `/usr/bin/orca` rejection.
- **Full Quantum Vibrational & Normal-Mode Parsers**: Parse harmonic frequencies, `km/mol` IR intensities, symmetry labels, imaginary mode counts, and `3N` Cartesian normal-mode displacements from **ORCA `.out`** (`VIBRATIONAL FREQUENCIES` + `IR SPECTRUM` + `NORMAL MODES`), **Gaussian `.log`** (including `HPModes` deduplication), **xTB `g98.out` / `vibspectrum`**, and **Molden `.molden`** files.
- **3D Normal-Mode Vibration Animation**: Play/pause live harmonic normal-mode vibrations in the 3D Viewer with adjustable amplitude and displacement vector arrows.
- **Experimental IR Overlay & Quantitative Matching**: Full IUPAC JCAMP-DX **ASDF** decoder (`AFFN`, `PAC`, `SQZ`, `DIF` checkpointing, `DUP`), `MICROMETERS` → `cm⁻¹` conversion, `%T` ↔ `Absorbance` preservation, 2-column `.csv`/`.tsv`/`.dpt`/`.txt` lab IR importer, optional NIST Chemistry WebBook `.jdx` lookup, live **HQI %** & **Pearson r** comparison, experimental peak picking, scroll-wheel zoom, and SVG/`.jdx`/`.csv` export.
- **Pure-Rust Explicit-H 3D Builder**: Automatically materializes 3D hydrogens + UFF pre-optimization for ORCA, Gaussian, and xTB input decks, plus a cached DFT deck preview and optional `--features ml-ir` backend trait seam.
- **MCP Server (19 Tools)**: Enriched `parse_quantum_output` with intensities & imaginary modes, and added `parse_jcamp_ir` and `run_xtb_ir`.

### Fixed
- **Core 2D Canvas Drawing UX**:
  - Ring (`Benzene`, `Cyclohexane`, etc.) and Chair tools now work on both click and drag (drag-to-rotate / drag-to-pick fusion side with anti-overlap flipping and bay atom auto-welding) with live hover and drag ghost previews.
  - Bond dragging auto-merges when the fixed-length/fixed-angle snapped endpoint lands on an existing atom, closing rings bond-by-bond cleanly with zero duplicate overlapping carbons.
  - Clicking atoms with Solid Bond now alternates trans zig-zag along chains and bisects the widest angular gap on 2-/3-/4-coordinate centers instead of curling or overlapping.
  - Acyclic Chain tool now attaches directly to existing atoms instead of stamping a disconnected fragment.
  - Clicking existing bonds with Wedge, HashWedge, HollowWedge, or Dative tools flips bond direction (`Begin ↔ End`); clicking a double bond with DoubleBond cycles alignment.
  - Double-clicking an atom with Marquee, Lasso, or Text opens inline atom editing immediately.
  - Dragging with Eraser sweeps and erases atoms/bonds in a single Undo step.
  - Added single-atom vertex dragging (with auto-weld on drop) and zoom-scaled screen-pixel hit tolerances (`atom_hit_tol` / `bond_hit_tol`).
- **Windows Release Packaging**: Excluded intermediate `.wixpdb` files from GitHub Release asset uploads.

## [1.0.0] - 2026-11-27

### Added
- **Publication-Grade 2D Drawing**: IUPAC 2008 (`GR-1`–`GR-10`) rendering, 16 clean-room journal stylesheets, 11 built-in template libraries (including Haworth sugars, macrocycles, bridged cages, lab glassware, BioGraphics, and organometallic catalysts), crossing-bond bridge clearance, smart reverse-order group typing, catalytic cycle layout, and live linked property/name captions.
- **Lossless Round-Trip Media (WordCraft / DeckCraft)**: Editable SVG (`<chiralcraft:document>`), Editable PNG (`tEXt` chunk), and Editable Vector PDF 1.4 with full extension persistence (energy profiles, table cell contents, TLC spot colors/labels/eluent), `--edit-svg` auto-save-on-close, and Wayland/X11 `image/svg+xml` clipboard support.
- **Interactive 3D & Crystal Viewer + UFF Force-Field**: Drag-to-orbit Ball & Stick / Stick / Space-Filling / Wireframe modes, `.xyz` / `.cif` / `.pdb` / VASP `POSCAR`/`CONTCAR` import with crystallographic unit-cell wireframe, pure-Rust UFF-lite 3D geometry optimization, click-to-measure distances/angles/dihedrals, and 3D rear-opacity depth cueing.
- **Computational Chemistry Bridge**: Export ORCA (`.inp`), Gaussian (`.gjf`), and xTB decks and parse `.out`/`.log` results (SCF energy, HOMO/LUMO gap, dipole moment, DFT IR frequencies, optimized 3D geometry).
- **Medicinal Chemistry, HRMS & Synthesis Suite**: Lipinski/Veber/Ghose/CNS MPO/QED/Fsp³/SA scorecard with live PAINS canvas highlights, exact 118-element HRMS isotope envelope & adduct simulator, physiological pH microspecies & pI, tautomer/stereoisomer/resonance enumerators, curved-arrow electron pushing, Markush R-group SAR grid, and offline retrosynthesis advisor.
- **Multi-Spectroscopy, Organometallics, Kinetics, ELN & LaTeX**: Simulated FT-IR (`%T`), UV–Vis $\lambda_{\max}$ & color swatch, 2D NMR (`COSY`/`HSQC`/`HMBC`), JCAMP-DX `.jdx` experimental NMR overlay, 18-electron Green CBC analyzer, Reaction Energy Profile (`ΔG‡`) & Eyring/Boltzmann kinetics, Polymer $M_n$/$M_w$/Đ & BigSMILES, ELN SI procedure generator + GHS safety handbook, `\chemfig`/TikZ LaTeX export, and multi-compound SDF/SMILES virtual screening.
- **Desktop Ergonomics & AI MCP Server**: Unified right-side tabbed Studio Dock, persistent `~/.config/chiralcraft/config.json`, `Ctrl+Scroll` zoom-to-cursor, middle-mouse pan, OS drag-and-drop, `:` / `F4` Command Palette, independent Light/Dark/Dark-Slide themes, Linux `.desktop` + MIME installer, and 17-tool stdio MCP server (`chiralcraft --mcp`).
