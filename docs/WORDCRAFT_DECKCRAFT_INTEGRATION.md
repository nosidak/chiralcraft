# WordCraft & DeckCraft Round-Trip Editable Graphic Integration

ChiralCraft natively embeds the full lossless chemical document inside exported **SVG**, **PNG**, and **PDF** files so office suites (**WordCraft**, **DeckCraft**, LibreOffice, etc.) can paste transparent publication graphics and double-click them later to edit in place.

---

## 1. Embedded Format Specification

### Editable SVG (`image/svg+xml`)
```xml
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:chiralcraft="https://chiralcraft.dev/ns/1.0"
     data-chiralcraft-version="1.0" ...>
  <metadata>
    <chiralcraft:document format="cdxml" encoding="base64">...</chiralcraft:document>
  </metadata>
  <!-- transparent bond-trimmed publication vector graphics -->
</svg>
```
- Transparent canvas background (zero white sticker circles around heteroatoms).
- Preserves molecules, stereochemistry, journal stylesheet, arrows, curves, tables (including cell text), TLC plates (including spot colors/labels/eluent), and reaction energy profiles (`ΔG‡`).

### Editable PNG (`image/png`)
- Contains a standard CRC32-checked PNG `tEXt` chunk with keyword `ChiralCraftDocument` and value `<base64(cdxml)>`.

### Editable Vector PDF (`application/pdf`)
- Pure PDF 1.4 vector paths & Helvetica text with trailing `%CHIRALCRAFT_CDXML_B64:<base64(cdxml)>`.

---

## 2. Clipboard Copy → Paste Workflow
1. User presses `Ctrl+C` (or `Edit > Copy As > Editable SVG`) in ChiralCraft.
2. ChiralCraft writes the editable SVG as `image/svg+xml` via `wl-copy` (Wayland) or `xclip` (X11).
3. WordCraft / DeckCraft receives `image/svg+xml` on `Ctrl+V` and stores the SVG inside the `.docx` / `.pptx` media archive unchanged.

---

## 3. Double-Click Round-Trip Editing Protocol (Host App Implementation)

When the user double-clicks an embedded SVG/PNG/PDF graphic in **WordCraft** or **DeckCraft**:

1. **Detect ChiralCraft payload** (fast substring check — zero external dependencies):
   - SVG: contains `<chiralcraft:document`
   - PNG: contains `ChiralCraftDocument\0`
   - PDF: contains `%CHIRALCRAFT_CDXML_B64:`
   - Or invoke `chiralcraft inspect-media <temp-file>` (exit code `0` and stdout `editable=true ...`).
2. **Write graphic bytes to a temp path**, e.g. `/tmp/wordcraft-chem-<uuid>.svg`.
3. **Launch ChiralCraft in round-trip mode**:
   ```bash
   chiralcraft --edit-svg /tmp/wordcraft-chem-<uuid>.svg
   ```
4. **Apply edits automatically**:
   - ChiralCraft displays a top banner (`Round-Trip Editable SVG`) with a **Save & Return** button.
   - Pressing `Ctrl+S` or closing the ChiralCraft window **atomically writes** the updated editable SVG back to `/tmp/wordcraft-chem-<uuid>.svg` (write temp + rename, so the host never sees a partial file).
   - Once the child process exits with status `0` (or via file watcher on save), WordCraft / DeckCraft reloads the file bytes and replaces the embedded image in place!

### Minimal Rust Host Helper (for WordCraft / DeckCraft)
```rust
use std::process::Command;

pub fn is_chiralcraft_svg(bytes: &[u8]) -> bool {
    bytes.windows(21).any(|w| w == b"<chiralcraft:document")
}

pub fn edit_chiralcraft_svg_blocking(svg_bytes: &[u8]) -> std::io::Result<Vec<u8>> {
    let tmp = std::env::temp_dir().join(format!("chiralcraft-rt-{}.svg", std::process::id()));
    std::fs::write(&tmp, svg_bytes)?;
    let status = Command::new("chiralcraft")
        .arg("--edit-svg")
        .arg(&tmp)
        .status()?;
    let updated = if status.success() { std::fs::read(&tmp)? } else { svg_bytes.to_vec() };
    let _ = std::fs::remove_file(&tmp);
    Ok(updated)
}
```
