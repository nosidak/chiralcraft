#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_SRC="$ROOT_DIR/target/release/chiralcraft"
BIN_DST="$HOME/.local/bin/chiralcraft"
APP_DIR="$HOME/.local/share/applications"
ICON_SVG_DIR="$HOME/.local/share/icons/hicolor/scalable/apps"
ICON_PNG_DIR="$HOME/.local/share/icons/hicolor/128x128/apps"
MIME_DIR="$HOME/.local/share/mime/packages"

if [[ ! -x "$BIN_SRC" ]]; then
  echo "Building ChiralCraft release binary..."
  (cd "$ROOT_DIR" && cargo build -p chiralcraft-app --release)
fi

mkdir -p "$HOME/.local/bin" "$APP_DIR" "$ICON_SVG_DIR" "$ICON_PNG_DIR" "$MIME_DIR"
cp -f "$BIN_SRC" "$BIN_DST"
chmod +x "$BIN_DST"

cp -f "$ROOT_DIR/assets/chiralcraft-icon.svg" "$ICON_SVG_DIR/chiralcraft.svg"
cp -f "$ROOT_DIR/assets/chiralcraft-icon-128.png" "$ICON_PNG_DIR/chiralcraft.png"
cp -f "$ROOT_DIR/packaging/chiralcraft-mime.xml" "$MIME_DIR/chiralcraft-mime.xml"

cat > "$APP_DIR/chiralcraft.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=ChiralCraft
GenericName=Chemical Structure & Cheminformatics Studio
Comment=Publication-grade 2D/3D chemical structure editor, spectroscopy, MedChem, and WordCraft/DeckCraft editable SVG studio
Exec=$BIN_DST %f
Icon=chiralcraft
Terminal=false
Categories=Science;Chemistry;Education;Graphics;
MimeType=chemical/x-cdxml;chemical/x-cdx;chemical/x-mdl-molfile;chemical/x-mdl-sdfile;chemical/x-xyz;chemical/x-cif;chemical/x-pdb;image/svg+xml;
StartupNotify=true
EOF

command -v update-mime-database >/dev/null 2>&1 && update-mime-database "$HOME/.local/share/mime" || true
command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database "$APP_DIR" || true
command -v gtk-update-icon-cache >/dev/null 2>&1 && gtk-update-icon-cache -f -t "$HOME/.local/share/icons/hicolor" 2>/dev/null || true

echo "Installed ChiralCraft to $BIN_DST"
echo "Desktop launcher: $APP_DIR/chiralcraft.desktop"
