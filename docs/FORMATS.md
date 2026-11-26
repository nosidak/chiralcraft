# File Formats & Lossless Round-Trip Matrix

ChiralCraft reads and writes **22+ chemical, crystallographic, spectral, quantum, LaTeX, and vector/raster formats** natively in pure Rust — with zero external converter binaries.

## Format Support Matrix

| Category | Format | Extensions | Read | Write | Lossless Scene Round-Trip | Notes |
|---|---|---|---|---|---|---|
| **Native Chemical** | CDXML | `.cdxml` | Yes | Yes | **Yes** | Preserves molecules, captions, arrows, curves, graphics, stylesheets, energy profiles, tables & TLC metadata |
| **Native Chemical** | Binary CDX | `.cdx`, `.cds` | Yes | Yes | **Yes** | Full binary object/property parser & serializer (22-byte & 28-byte headers) |
| **Lossless Media** | Editable SVG | `.svg` | Yes | Yes | **Yes** | Embeds base64 CDXML inside `<metadata><chiralcraft:document>` — transparent & dark-slide modes |
| **Lossless Media** | Editable PNG | `.png` | Yes | Yes | **Yes** | Configurable DPI (`150/300/600`), embeds lossless CDXML inside PNG `tEXt` chunk `ChiralCraftDocument` |
| **Lossless Media** | Editable Vector PDF | `.pdf` | Yes | Yes | **Yes** | Native pure-Rust PDF 1.4 vector graphics with embedded `%CHIRALCRAFT_CDXML_B64:` payload |
| **MDL / CTAB** | Molfile V2000 & V3000 | `.mol` | Yes | Yes | — | Stereo wedges, charges, isotopes, radicals, Data SGroups (`DAT`), Enhanced Stereo (`CFG`/`STEABS`/`STEREAC`) |
| **MDL / CTAB** | Structure-Data File | `.sdf`, `.sd` | Yes | Yes | — | Multi-record libraries, property tags, virtual screening & Markush export |
| **MDL / CTAB** | Reaction File | `.rxn` | Yes | Yes | — | Reactants, products, atom-map numbers & reaction centers (V2000/V3000) |
| **Line Notations** | SMILES | `.smi`, `.smiles`, CLI | Yes | Yes | — | Aromaticity, stereochemistry (`@`/`@@`, `/`/`\`), isotopes, charges, disconnected salts |
| **Line Notations** | InChI & InChIKey | `.inchi` | Dialog | Yes | — | Standard InChI + SHA-256 InChIKey generation & dialog paste |
| **Line Notations** | SLN | `.sln` | — | Yes | — | SYBYL Line Notation graph traversal serializer |
| **Interchange** | Chemical Markup Language | `.cml` | Yes | Yes | — | XML atom/bond array interchange |
| **Interchange** | Connection Table | `.ct` | Yes | Yes | — | Coordinate + connectivity table with wedge/hash codes |
| **Interchange** | Trivial Graph / Sketch | `.tgf`, `.skc` | Yes | Yes | — | Legacy chemical graph & binary sketch interchange |
| **3D & Crystals** | XYZ Coordinates | `.xyz` | Yes | Yes* | — | Cartesian 3D coordinates + automatic covalent bond perception (*exported via xTB deck) |
| **3D & Crystals** | Crystallographic Info | `.cif` | Yes | — | — | Fractional coordinates, unit cell (`a, b, c, α, β, γ`), live unit-cell wireframe in 3D Viewer |
| **3D & Crystals** | Protein Data Bank | `.pdb` | Yes | — | — | `ATOM`/`HETATM` & `CONECT` records with 3D coordinates |
| **3D & Crystals** | VASP Crystal | `POSCAR`, `CONTCAR`, `.vasp` | Yes | — | — | Direct & Cartesian lattice vectors + periodic unit-cell wireframe |
| **Quantum** | ORCA / Gaussian / xTB | `.inp`, `.gjf`, `.out`, `.log` | Yes | Yes | — | Export Opt+Freq / SP / TD-DFT / GIAO-NMR decks; parse SCF energy, HOMO–LUMO, dipole, IR & 3D coords |
| **Spectroscopy** | JCAMP-DX | `.jdx`, `.dx` | Yes | — | — | Experimental NMR spectrum overlay (`XYDATA` & `PEAK TABLE`) |
| **Biopolymers** | HELM & FASTA | Clipboard / Dialog | Yes | Yes | — | Peptide, DNA, and RNA biopolymer notation |
| **Publishing & 3D** | LaTeX / TikZ | `.tex` | — | Yes | — | Compilable `\chemfig{...}`, `\schemestart...\schemestop`, and `\begin{tikzpicture}` |
| **Publishing & 3D** | Vector Legacy & 3D Print | `.eps`, `.emf`, `.3mf` | — | Yes | — | Encapsulated PostScript, Windows Enhanced Metafile, and 3MF Ball-and-Stick / Stick meshes |

## Extension Persistence inside CDXML, SVG, PNG & PDF

Beyond standard CDXML molecules, arrows, captions, curves, and graphics, ChiralCraft serializes rich studio extensions inside a lossless `<objecttag Name="chiralcraft_extensions">` payload so that reopening an `.svg`, `.png`, `.pdf`, `.cdxml`, or `.cdx` file restores:
- **Reaction Energy Profiles (`ΔG‡`)**: Pathway A & B states, energies, transition-state markers, temperature, units (`kcal/mol` or `kJ/mol`), and canvas bounding box.
- **Data Tables**: Full cell text grid, row/column dimensions, and placement.
- **TLC Plates**: Eluent label, spot `Rf` values, lane labels, and custom RGB spot colors.
