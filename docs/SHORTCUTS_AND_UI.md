# Keyboard Shortcuts, Context Hotkeys & UI Guide

ChiralCraft is designed for fast, fluid, keyboard-assisted structure drawing alongside full menu and palette control.

---

## 1. Command Palette (`:` or `F4`)

Press `:` (when not editing text) or `F4` anywhere to open the fuzzy Command Palette. In addition to fuzzy-searching every action, you can type direct parameterized commands:

| Command | Action |
|---|---|
| `rotate 30` / `rot -45` | Rotate selected (or all) structures by any angle in degrees |
| `scale 150` | Scale selected structures by percentage |
| `smiles CC(=O)Oc1ccccc1C(=O)O` | Parse SMILES, clean up 2D layout, and stamp onto canvas |
| `name caffeine` | Resolve IUPAC/common name and stamp onto canvas |
| `style nature` / `style acs` / `style rsc` | Apply a matching journal stylesheet immediately |
| `theme light` / `theme dark` / `theme slide` | Switch UI and canvas theme (`ClassicLight`, `DarkChrome`, `DarkSlide`) |
| `nmr stamp` / `nmr 13c` / `nmr clear` | Stamp ¹H or ¹³C NMR shift badges on canvas, or clear them |
| `cycle circle` / `cycle ellipse` / `cycle starter` | Arrange selected molecules as a catalytic cycle or insert starter 4-step cycle |
| `core adamantane` / `core haworth` / `core porphine` / `core crown` | Stamp special chemical core directly |
| `dpi 600` | Update editable PNG export resolution (`150` / `300` / `600`) |
| `3d` · `uff` · `diff` · `medchem` · `hrms` · `retro` · `ir` · `18e` · `energy` · `latex` · `quantum` · `dock` | Open any studio or trigger instant UFF 3D optimization |

---

## 2. Global Menu & Canvas Shortcuts

| Shortcut | Action |
|---|---|
| `Ctrl+N` / `Ctrl+O` / `Ctrl+S` / `Shift+Ctrl+S` / `Ctrl+W` | New · Open · Save · Save As · Close tab |
| `Ctrl+Z` / `Shift+Ctrl+Z` / `Ctrl+Y` | Undo · Redo · Repeat last action |
| `Ctrl+X` / `Ctrl+C` / `Ctrl+V` / `Del` | Cut · Copy (also pushes Editable SVG to OS clipboard) · Paste · Clear |
| `Ctrl+A` / `Shift+Ctrl+I` | Select All · Invert Selection |
| `Shift+Ctrl+K` / `Shift+Ctrl+X` | Clean Up Structure · Clean Up Reaction |
| `Alt+K` | Toggle Aromatic Ring Circle Display |
| `Shift+Ctrl+N` / `Alt+Ctrl+N` | Convert Name to Structure · Convert Structure to Name |
| `Shift+Ctrl+F` | Search PubChem & ChEMBL Database |
| `Alt+Ctrl+C` / `Alt+Ctrl+P` | Copy as SMILES · Paste SMILES |
| `Ctrl+L` / `Ctrl+E` | Toggle Fixed Lengths · Toggle Fixed Angles (`15°` steps) |
| `Ctrl+G` / `Shift+Ctrl+G` | Group · Ungroup selected structures |
| `F2` / `F3` | Bring to Front · Send to Back |
| `F5` / `F6` / `F7` / `F8` | Actual Size (100%) · Fit Document · Magnify · Reduce |
| `Ctrl+Wheel` / Middle-drag | Zoom to cursor · Pan canvas |
| `F11` / `Ctrl+H` / `F1` | Toggle Rulers · Toggle Crosshair · Shortcuts Cheat Sheet |

---

## 3. Context-Sensitive Atom & Bond Hotkeys (98 Bindings)

Hover any **atom** or **bond** on the canvas and press a single key — no tool switching needed:

### Hovering an Atom
- **Elements & Groups**: `c` (C), `C` (Cl), `n` (N), `o` (O), `s` (S), `p` (P), `f` (F), `F` (CF₃), `b` (B), `B` (Br), `i` (I), `h` (H), `k` (K), `l` (Li), `m` (Me), `a` (Ac), `t` (tBu), `P` (Ph), `S` (Si), `M` (Mg), `Z` (Zn)
- **Quick Sprouting**: `1` (single C–C sprout), `2` (`=O` carbonyl sprout), `3` (≡N nitrile / alkyne sprout), `4` (4-ring), `5` (cyclopentane), `6` (cyclohexane), `a` / `r` (aromatic/ring actions)
- **Charges & Editing**: `+` (charge `+1`), `-` (charge `-1`), `Enter` (open inline atom label editor), `Backspace` / `Delete` (erase atom)

### Smart Reverse Functional-Group Typing
Press `Enter` on any atom to type a custom label or nickname. Left-facing reversed groups are automatically recognized and remain expandable via `Structure › Expand Label`:
`HO` · `H2N` · `HS` · `O2N` (`O₂N`) · `MeO` · `EtO` · `iPrO` · `tBuO` · `AcO` · `PhO` · `BnO` · `TBSO` · `TMSO` · `TsO` · `MsO` · `TfO` · `F3C` (`F₃C`) · `Cl3C` · `NC` · `HO2C` (`HO₂C`) · `MeO2C` · `EtO2C` · `BocHN` · `CbzHN` · `FmocHN` · `AcHN` · `MeHN` · `Me2N` · `Et2N` · `O2S`

### Hovering a Bond
- `1` / `2` / `3` — Set Single / Double / Triple bond
- `w` / `h` / `b` / `d` / `y` — Wedge / Hashed-Wedge / Bold / Dashed / Wavy stereobond
- `4` / `5` / `6` / `7` / `8` — Fuse a 4-, 5-, 6-, 7-, or 8-membered ring directly across the hovered bond
- `f` — Flip stereobond wedge direction; `Backspace` / `Delete` — Erase bond

---

## 4. Drawing Interactions (Mouse)

- **Click empty canvas** with any Bond, Ring, Arrow, Shape, Bracket, Curve, Chain, Orbital, Symbol, Table, TLC, or Template tool to stamp a default-sized object immediately.
- **Click an existing atom** with a Ring tool to perform **Spiro Ring Fusion**; **click a bond midpoint** to perform **Edge Ring Fusion** (e.g., Benzene + Benzene → Naphthalene).
- **Drag a bond** from one atom onto an atom in a different structure to automatically **merge** the two molecules with a single bond.
- **Right-click** any atom or bond for quick properties, label expansion, stereobond flipping, and deletion.
