#!/usr/bin/env bash
set -euo pipefail

VERSION="${1:?Usage: package.sh <version> <arm64_bin> <x64_bin> [out_dir]}"
ARM64_BIN="${2:?Missing arm64 binary path}"
X64_BIN="${3:?Missing x86_64 binary path}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${4:-$ROOT_DIR/dist}"

mkdir -p "$OUT_DIR"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

# 1. Per-arch standalone tarballs
for PAIR in "aarch64:$ARM64_BIN" "x86_64:$X64_BIN"; do
  ARCH="${PAIR%%:*}"
  BIN="${PAIR#*:}"
  NAME="chiralcraft-${VERSION}-macos-${ARCH}"
  STAGE="$WORK_DIR/$NAME"
  mkdir -p "$STAGE"
  cp -f "$BIN" "$STAGE/chiralcraft"
  chmod +x "$STAGE/chiralcraft"
  cp -f "$ROOT_DIR/README.md" "$ROOT_DIR/LICENSE" "$ROOT_DIR/CHANGELOG.md" "$STAGE/"
  tar -C "$WORK_DIR" -czf "$OUT_DIR/${NAME}.tar.gz" "$NAME"
done

# 2. Universal fat binary (arm64 + x86_64)
UNIV_BIN="$WORK_DIR/chiralcraft-universal"
lipo -create -output "$UNIV_BIN" "$ARM64_BIN" "$X64_BIN"
chmod +x "$UNIV_BIN"
lipo -info "$UNIV_BIN"

# 3. Generate macOS .icns icon from PNG
ICONSET="$WORK_DIR/ChiralCraft.iconset"
mkdir -p "$ICONSET"
SRC_PNG="$ROOT_DIR/assets/chiralcraft-icon.png"
for SZ in 16 32 128 256 512; do
  sips -z "$SZ" "$SZ" "$SRC_PNG" --out "$ICONSET/icon_${SZ}x${SZ}.png" >/dev/null
  SZ2=$((SZ * 2))
  sips -z "$SZ2" "$SZ2" "$SRC_PNG" --out "$ICONSET/icon_${SZ}x${SZ}@2x.png" >/dev/null
done
iconutil -c icns "$ICONSET" -o "$WORK_DIR/ChiralCraft.icns"

# 4. Assemble Universal ChiralCraft.app bundle
APP="$WORK_DIR/ChiralCraft.app"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp -f "$UNIV_BIN" "$APP/Contents/MacOS/chiralcraft"
cp -f "$WORK_DIR/ChiralCraft.icns" "$APP/Contents/Resources/ChiralCraft.icns"
cat > "$APP/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleDevelopmentRegion</key>
  <string>en</string>
  <key>CFBundleDisplayName</key>
  <string>ChiralCraft</string>
  <key>CFBundleExecutable</key>
  <string>chiralcraft</string>
  <key>CFBundleIconFile</key>
  <string>ChiralCraft</string>
  <key>CFBundleIdentifier</key>
  <string>io.github.nosidak.chiralcraft</string>
  <key>CFBundleInfoDictionaryVersion</key>
  <string>6.0</string>
  <key>CFBundleName</key>
  <string>ChiralCraft</string>
  <key>CFBundlePackageType</key>
  <string>APPL</string>
  <key>CFBundleShortVersionString</key>
  <string>${VERSION}</string>
  <key>CFBundleVersion</key>
  <string>${VERSION}</string>
  <key>LSApplicationCategoryType</key>
  <string>public.app-category.education</string>
  <key>LSMinimumSystemVersion</key>
  <string>11.0</string>
  <key>NSHighResolutionCapable</key>
  <true/>
  <key>CFBundleDocumentTypes</key>
  <array>
    <dict>
      <key>CFBundleTypeName</key>
      <string>Chemical Structure Document</string>
      <key>CFBundleTypeRole</key>
      <string>Editor</string>
      <key>CFBundleTypeExtensions</key>
      <array>
        <string>cdxml</string>
        <string>cdx</string>
        <string>mol</string>
        <string>sdf</string>
        <string>xyz</string>
        <string>cif</string>
        <string>pdb</string>
        <string>svg</string>
      </array>
    </dict>
  </array>
</dict>
</plist>
EOF

# Ad-hoc sign so Apple Silicon Gatekeeper allows local execution cleanly
codesign --force --deep --sign - "$APP"

# 5. Universal .app .zip
(cd "$WORK_DIR" && ditto -c -k --sequesterRsrc --keepParent "ChiralCraft.app" \
  "$OUT_DIR/chiralcraft-${VERSION}-macos-universal.zip")

# 6. Universal drag-to-Applications .dmg
DMG_STAGE="$WORK_DIR/dmg-stage"
mkdir -p "$DMG_STAGE"
cp -R "$APP" "$DMG_STAGE/"
ln -s /Applications "$DMG_STAGE/Applications"
cp -f "$ROOT_DIR/README.md" "$ROOT_DIR/LICENSE" "$DMG_STAGE/"
hdiutil create \
  -volname "ChiralCraft ${VERSION}" \
  -srcfolder "$DMG_STAGE" \
  -ov -format UDZO \
  "$OUT_DIR/chiralcraft-${VERSION}-macos-universal.dmg"

echo "==> macOS packaging complete:"
ls -lh "$OUT_DIR/chiralcraft-${VERSION}-macos-"*
