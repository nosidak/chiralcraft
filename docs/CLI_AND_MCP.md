# CLI & MCP Server Reference

ChiralCraft bundles a complete headless scientific CLI and a zero-dependency stdio **Model Context Protocol (MCP)** server directly into the single `chiralcraft` binary.

---

## 1. Command-Line Interface (CLI)

Every chemistry subcommand accepts **either a file path** (any supported format) **or a raw SMILES string**.

```bash
chiralcraft --help                      # Full usage help
chiralcraft --version                   # Print version (chiralcraft 1.1.1)
chiralcraft [FILE]                      # Launch GUI (optionally opening FILE)
chiralcraft --edit-svg <FILE.svg>       # Round-trip edit an embedded ChiralCraft SVG
chiralcraft --screenshot <PRESET> <OUT> # Capture deterministic GPU screenshot (hero, medchem, spectroscopy, 3d-quantum, synthesis, organometallic, dark-slide)
chiralcraft --mcp                       # Start stdio JSON-RPC 2.0 MCP server
```

### Scientific & Conversion Subcommands

| Command | Usage | Output |
|---|---|---|
| `analyze` | `chiralcraft analyze <FILE\|SMILES>` | Formula, Exact Mass, MW, Elemental Analysis, BP, MP, LogP, MR, tPSA, LogS, pKa |
| `convert` | `chiralcraft convert <IN> <OUT>` | Convert between any readable input (or SMILES) and any writable output format |
| `diff` | `chiralcraft diff <LEFT> <RIGHT>` | Chemistry-aware diff across `.cdxml`, `.cdx`, `.svg`, `.png`, `.pdf`, `.mol`, `.sdf` |
| `inspect-svg` / `inspect-media` | `chiralcraft inspect-media <FILE>` | Verifies embedded lossless ChiralCraft payload in `.svg`, `.png`, or `.pdf` (exits `0` if editable, `2` otherwise) |
| `nmr` | `chiralcraft nmr <FILE\|SMILES> [--13c]` | Predicted ¹H-NMR (400 MHz) or ¹³C-NMR (100 MHz) chemical shifts, multiplicities & proton counts |
| `ir` | `chiralcraft ir <FILE\|SMILES>` | Predicted FT-IR bands (`cm⁻¹`, strength, functional assignment) + UV–Vis `λmax` & color |
| `medchem` | `chiralcraft medchem <FILE\|SMILES>` | Lipinski, Veber, Ghose, Pfizer CNS MPO, QED, Fsp³, SA score & PAINS/Brenk alerts |
| `hrms` | `chiralcraft hrms <FILE\|SMILES> [--adduct <A>]` | Neutral exact mass, exact isotope peaks (`m/z`, rel%) for selected adduct (default `[M+H]+`) |
| `retro` | `chiralcraft retro <FILE\|SMILES>` | Ranked retrosynthetic disconnections, synthons & forward reagents |
| `18e` | `chiralcraft 18e <FILE\|SMILES>` | Transition-metal 18-electron & Green CBC (`MLxXy`) oxidation state, `dⁿ`, geometry & status |
| `export-tex` | `chiralcraft export-tex <FILE\|SMILES>` | Compilable LaTeX `\schemestart` / `\chemfig` snippet |
| `screen` | `chiralcraft screen <LIB.smi\|csv\|sdf>` | Batch MedChem, CNS MPO, QED & PAINS screening table |
| `quantum-in` | `chiralcraft quantum-in <FILE\|SMILES> [--engine orca\|gaussian\|xtb]` | Lifts 2D→3D and prints a ready-to-run ORCA, Gaussian, or GFN2-xTB input deck |
| `quantum-parse` | `chiralcraft quantum-parse <JOB.out\|log>` | Extracts convergence, SCF energy (`Eh`), HOMO–LUMO gap (`eV`), dipole (`D`) & IR frequency count |

---

## 2. Built-In MCP Server (`chiralcraft --mcp`)

ChiralCraft implements **Model Context Protocol (JSON-RPC 2.0 over newline-delimited stdio)** with **19 chemistry and publishing tools** so AI assistants and coding agents can draw, verify, analyze, and render chemistry directly.

### Client Configuration Example

Add ChiralCraft to any MCP-compatible host (JSON config):

```json
{
  "mcpServers": {
    "chiralcraft": {
      "command": "chiralcraft",
      "args": ["--mcp"]
    }
  }
}
```

### All 19 Exposed MCP Tools

| # | Tool Name | Required Args | Optional Args | Description |
|---|---|---|---|---|
| 1 | `analyze_smiles` | `smiles` | — | Formula, exact mass, MW, elemental analysis, LogP, tPSA, LogS, pKa, IUPAC name & warnings |
| 2 | `predict_nmr` | `smiles` | `nucleus` (`"1H"` \| `"13C"`) | Predict ¹H-NMR (400 MHz) or ¹³C-NMR (100 MHz) shifts and multiplicities |
| 3 | `name_to_structure` | `name` | `stylesheet` | Convert an IUPAC/common name to SMILES, MOL V3000, and editable ChiralCraft SVG |
| 4 | `render_editable_svg` | `smiles` | `stylesheet`, `transparent`, `dark_ink`, `caption`, `output_path` | Render SMILES to a self-contained lossless Editable SVG (optionally written to disk) |
| 5 | `search_chem_db` | `query` | `source` (`"both"` \| `"pubchem"` \| `"chembl"`) | Live query against PubChem and/or ChEMBL by name, CAS, CID, ChEMBL ID, formula, or SMILES |
| 6 | `inspect_file` | `path` | — | Inspect any chemical or editable media file (`.cdxml`, `.cdx`, `.mol`, `.sdf`, `.xyz`, `.cif`, `.pdb`, `POSCAR`, `.svg`, `.png`, `.pdf`) |
| 7 | `convert_file` | `input_path`, `output_path` | — | Convert between supported chemical, 3D, and vector/raster formats |
| 8 | `diff_files` | `left_path`, `right_path` | — | Perform a chemistry-aware structural & stereochemical diff between two files |
| 9 | `medchem_scorecard` | `smiles` | — | Compute Lipinski, Veber, Ghose, CNS MPO, QED, Fsp³, SA score & PAINS/Brenk alerts |
| 10 | `simulate_hrms` | `smiles` | `adduct` | Simulate exact HRMS isotopic envelope peaks and ESI+/ESI− adduct table |
| 11 | `retrosynthesis_suggest` | `smiles` | — | Suggest ranked strategic retrosynthetic disconnections, synthons & reagents |
| 12 | `predict_spectroscopy` | `smiles` | — | Predict FT-IR bands, UV–Vis `λmax`/color, and 2D NMR (COSY/HSQC/HMBC) cross-peak counts |
| 13 | `analyze_organometallic` | `smiles` | — | Analyze transition-metal complexes via the 18-electron rule & Green CBC (`MLxXy`) |
| 14 | `export_latex` | `smiles` | — | Export compilable LaTeX `\chemfig{...}` and `tikzpicture` vector code |
| 15 | `generate_si_procedure` | `smiles` | `name`, `solvent`, `yield_pct` | Generate ACS/JOC Supporting Information experimental procedure + NMR/IR/HRMS block |
| 16 | `export_quantum_input` | `smiles` | `program` (`"orca"` \| `"gaussian"` \| `"xtb"`), `method`, `basis`, `solvent` | Generate ready-to-run quantum chemistry input decks from SMILES |
| 17 | `parse_quantum_output` | `path` | — | Parse ORCA `.out`, Gaussian `.log`, xTB `g98.out`/`vibspectrum`, or Molden files for SCF energy, HOMO–LUMO gap, dipole & IR modes |
| 18 | `parse_jcamp_ir` | `path` | `compare_smiles` | Parse experimental FT-IR JCAMP-DX (`.jdx`/`.dx`, ASDF/SQZ/DIF/DUP) or CSV/DPT files, pick peaks & compute HQI/Pearson match |
| 19 | `run_xtb_ir` | `smiles` | `method` (`"gfn2"` \| `"gfnff"`), `solvent` | Run local GFN2-xTB or GFN-FF harmonic IR (`--ohess`) and return scaled wavenumbers & intensities (km/mol) |
