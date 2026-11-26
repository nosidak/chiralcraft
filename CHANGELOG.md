# Changelog

All notable changes to **ChiralCraft** are documented in this file.

## [1.0.0] - 2026-11-27

### Added
- **Publication-Grade 2D Drawing**: IUPAC 2008 (`GR-1`–`GR-10`) rendering, 16 clean-room journal stylesheets, 11 built-in template libraries (including Haworth sugars, macrocycles, bridged cages, lab glassware, BioGraphics, and organometallic catalysts), crossing-bond bridge clearance, smart reverse-order group typing, catalytic cycle layout, and live linked property/name captions.
- **Lossless Round-Trip Media (WordCraft / DeckCraft)**: Editable SVG (`<chiralcraft:document>`), Editable PNG (`tEXt` chunk), and Editable Vector PDF 1.4 with full extension persistence (energy profiles, table cell contents, TLC spot colors/labels/eluent), `--edit-svg` auto-save-on-close, and Wayland/X11 `image/svg+xml` clipboard support.
- **Interactive 3D & Crystal Viewer + UFF Force-Field**: Drag-to-orbit Ball & Stick / Stick / Space-Filling / Wireframe modes, `.xyz` / `.cif` / `.pdb` / VASP `POSCAR`/`CONTCAR` import with crystallographic unit-cell wireframe, pure-Rust UFF-lite 3D geometry optimization, click-to-measure distances/angles/dihedrals, and 3D rear-opacity depth cueing.
- **Computational Chemistry Bridge**: Export ORCA (`.inp`), Gaussian (`.gjf`), and xTB decks and parse `.out`/`.log` results (SCF energy, HOMO/LUMO gap, dipole moment, DFT IR frequencies, optimized 3D geometry).
- **Medicinal Chemistry, HRMS & Synthesis Suite**: Lipinski/Veber/Ghose/CNS MPO/QED/Fsp³/SA scorecard with live PAINS canvas highlights, exact 118-element HRMS isotope envelope & adduct simulator, physiological pH microspecies & pI, tautomer/stereoisomer/resonance enumerators, curved-arrow electron pushing, Markush R-group SAR grid, and offline retrosynthesis advisor.
- **Multi-Spectroscopy, Organometallics, Kinetics, ELN & LaTeX**: Simulated FT-IR (`%T`), UV–Vis $\lambda_{\max}$ & color swatch, 2D NMR (`COSY`/`HSQC`/`HMBC`), JCAMP-DX `.jdx` experimental NMR overlay, 18-electron Green CBC analyzer, Reaction Energy Profile (`ΔG‡`) & Eyring/Boltzmann kinetics, Polymer $M_n$/$M_w$/Đ & BigSMILES, ELN SI procedure generator + GHS safety handbook, `\chemfig`/TikZ LaTeX export, and multi-compound SDF/SMILES virtual screening.
- **Desktop Ergonomics & AI MCP Server**: Unified right-side tabbed Studio Dock, persistent `~/.config/chiralcraft/config.json`, `Ctrl+Scroll` zoom-to-cursor, middle-mouse pan, OS drag-and-drop, `:` / `F4` Command Palette, independent Light/Dark/Dark-Slide themes, Linux `.desktop` + MIME installer, and 17-tool stdio MCP server (`chiralcraft --mcp`).
