#!/usr/bin/env bash
set -euo pipefail

VERSION="${1:?Usage: package.sh <version> <arch: x86_64|aarch64> [bin_path] [out_dir]}"
ARCH="${2:?Usage: package.sh <version> <arch: x86_64|aarch64> [bin_path] [out_dir]}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN_PATH="${3:-$ROOT_DIR/target/release/chiralcraft}"
OUT_DIR="${4:-$ROOT_DIR/dist}"

if [[ ! -x "$BIN_PATH" ]]; then
  echo "error: binary not found or not executable: $BIN_PATH" >&2
  exit 1
fi

case "$ARCH" in
  x86_64)
    DEB_ARCH="amd64"
    RPM_ARCH="x86_64"
    PAC_ARCH="x86_64"
    APPIMAGE_ARCH="x86_64"
    ;;
  aarch64)
    DEB_ARCH="arm64"
    RPM_ARCH="aarch64"
    PAC_ARCH="aarch64"
    APPIMAGE_ARCH="aarch64"
    ;;
  *)
    echo "error: unsupported arch '$ARCH' (expected x86_64 or aarch64)" >&2
    exit 1
    ;;
esac

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PKG_BASE="chiralcraft-${VERSION}-linux-${ARCH}"

# Shared FHS staging tree
STAGE="$WORK_DIR/stage"
install -Dm755 "$BIN_PATH" "$STAGE/usr/bin/chiralcraft"
install -Dm644 "$ROOT_DIR/assets/chiralcraft.desktop" "$STAGE/usr/share/applications/chiralcraft.desktop"
install -Dm644 "$ROOT_DIR/assets/chiralcraft-icon-128.png" "$STAGE/usr/share/icons/hicolor/128x128/apps/chiralcraft.png"
install -Dm644 "$ROOT_DIR/assets/chiralcraft-icon.svg" "$STAGE/usr/share/icons/hicolor/scalable/apps/chiralcraft.svg"
install -Dm644 "$ROOT_DIR/packaging/chiralcraft-mime.xml" "$STAGE/usr/share/mime/packages/chiralcraft-mime.xml"
install -Dm644 "$ROOT_DIR/packaging/linux/io.github.nosidak.chiralcraft.metainfo.xml" \
  "$STAGE/usr/share/metainfo/io.github.nosidak.chiralcraft.metainfo.xml"
install -Dm644 "$ROOT_DIR/LICENSE" "$STAGE/usr/share/licenses/chiralcraft/LICENSE"
install -Dm644 "$ROOT_DIR/README.md" "$STAGE/usr/share/doc/chiralcraft/README.md"

# 1. Portable .tar.gz
echo "==> Building ${PKG_BASE}.tar.gz"
TAR_ROOT="$WORK_DIR/${PKG_BASE}"
mkdir -p "$TAR_ROOT/assets" "$TAR_ROOT/packaging" "$TAR_ROOT/scripts"
install -m755 "$BIN_PATH" "$TAR_ROOT/chiralcraft"
cp -f "$ROOT_DIR/README.md" "$ROOT_DIR/LICENSE" "$ROOT_DIR/CHANGELOG.md" "$TAR_ROOT/"
cp -f "$ROOT_DIR/assets/chiralcraft.desktop" "$ROOT_DIR/assets/chiralcraft-icon-128.png" "$ROOT_DIR/assets/chiralcraft-icon.svg" "$TAR_ROOT/assets/"
cp -f "$ROOT_DIR/packaging/chiralcraft-mime.xml" "$TAR_ROOT/packaging/"
cp -f "$ROOT_DIR/scripts/install-linux-desktop.sh" "$TAR_ROOT/scripts/"
tar -C "$WORK_DIR" -czf "$OUT_DIR/${PKG_BASE}.tar.gz" "${PKG_BASE}"

# 2. Debian / Ubuntu .deb
if command -v dpkg-deb >/dev/null 2>&1; then
  echo "==> Building ${PKG_BASE}.deb"
  DEB_ROOT="$WORK_DIR/deb"
  cp -a "$STAGE" "$DEB_ROOT"
  INSTALLED_KB="$(du -sk "$DEB_ROOT/usr" | awk '{print $1}')"
  mkdir -p "$DEB_ROOT/DEBIAN"
  cat > "$DEB_ROOT/DEBIAN/control" <<EOF
Package: chiralcraft
Version: ${VERSION}
Section: science
Priority: optional
Architecture: ${DEB_ARCH}
Installed-Size: ${INSTALLED_KB}
Maintainer: nosidak <https://github.com/nosidak/chiralcraft>
Homepage: https://github.com/nosidak/chiralcraft
Depends: libc6, libgcc-s1, libgl1 | libegl1, libxkbcommon0
Description: Pure-Rust chemical structure editor and cheminformatics studio
 ChiralCraft is an open-source 2D/3D chemical structure editor, spectroscopy
 simulator, MedChem ADMET/PAINS workbench, retrosynthesis advisor, CLI batch
 converter, and MCP server licensed under GNU AGPL-3.0-or-later.
EOF
  cat > "$DEB_ROOT/DEBIAN/postinst" <<'EOF'
#!/bin/sh
set -e
command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database -q /usr/share/applications || true
command -v update-mime-database >/dev/null 2>&1 && update-mime-database /usr/share/mime || true
command -v gtk-update-icon-cache >/dev/null 2>&1 && gtk-update-icon-cache -q -f -t /usr/share/icons/hicolor || true
EOF
  chmod 755 "$DEB_ROOT/DEBIAN/postinst"
  dpkg-deb --root-owner-group --build "$DEB_ROOT" "$OUT_DIR/${PKG_BASE}.deb"
fi

# 3. Fedora / RHEL / openSUSE .rpm
if command -v rpmbuild >/dev/null 2>&1; then
  echo "==> Building ${PKG_BASE}.rpm"
  RPM_TOP="$WORK_DIR/rpmbuild"
  mkdir -p "$RPM_TOP"/{BUILD,RPMS,SOURCES,SPECS,SRPMS}
  cat > "$RPM_TOP/SPECS/chiralcraft.spec" <<EOF
%global debug_package %{nil}
%global _build_id_links none

Name:           chiralcraft
Version:        ${VERSION}
Release:        1%{?dist}
Summary:        Pure-Rust chemical structure editor and cheminformatics studio
License:        AGPL-3.0-or-later
URL:            https://github.com/nosidak/chiralcraft
BuildArch:      ${RPM_ARCH}

%description
ChiralCraft is a free, open-source, pure-Rust desktop chemical structure editor,
spectroscopy simulator, MedChem workbench, 3D viewer, CLI converter, and MCP server.

%install
cp -a ${STAGE}/* %{buildroot}/

%files
/usr/bin/chiralcraft
/usr/share/applications/chiralcraft.desktop
/usr/share/icons/hicolor/128x128/apps/chiralcraft.png
/usr/share/icons/hicolor/scalable/apps/chiralcraft.svg
/usr/share/mime/packages/chiralcraft-mime.xml
/usr/share/metainfo/io.github.nosidak.chiralcraft.metainfo.xml
/usr/share/licenses/chiralcraft/LICENSE
/usr/share/doc/chiralcraft/README.md
EOF
  rpmbuild --define "_topdir $RPM_TOP" --target "$RPM_ARCH" -bb "$RPM_TOP/SPECS/chiralcraft.spec"
  BUILT_RPM="$(find "$RPM_TOP/RPMS" -name '*.rpm' | head -n 1)"
  cp -f "$BUILT_RPM" "$OUT_DIR/${PKG_BASE}.rpm"
fi

# 4. Arch Linux .pkg.tar.zst
if command -v zstd >/dev/null 2>&1; then
  echo "==> Building ${PKG_BASE}.pkg.tar.zst"
  PAC_ROOT="$WORK_DIR/pacman"
  cp -a "$STAGE" "$PAC_ROOT"
  PKG_BYTES="$(du -sb "$PAC_ROOT/usr" | awk '{print $1}')"
  BUILD_DATE="$(date +%s)"
  cat > "$PAC_ROOT/.PKGINFO" <<EOF
pkgname = chiralcraft
pkgbase = chiralcraft
pkgver = ${VERSION}-1
pkgdesc = Pure-Rust chemical structure editor and cheminformatics studio
url = https://github.com/nosidak/chiralcraft
builddate = ${BUILD_DATE}
packager = nosidak <https://github.com/nosidak/chiralcraft>
size = ${PKG_BYTES}
arch = ${PAC_ARCH}
license = AGPL-3.0-or-later
depend = glibc
depend = gcc-libs
depend = libxkbcommon
depend = libglvnd
EOF
  (cd "$PAC_ROOT" && tar --owner=0 --group=0 -cf - .PKGINFO usr | zstd -19 -q -o "$OUT_DIR/${PKG_BASE}.pkg.tar.zst")
fi

# 5. Portable .AppImage (when appimagetool is available on PATH or APPIMAGETOOL env var)
APPIMAGETOOL="${APPIMAGETOOL:-$(command -v appimagetool || true)}"
if [[ -n "$APPIMAGETOOL" && -x "$APPIMAGETOOL" ]]; then
  echo "==> Building ${PKG_BASE}.AppImage"
  APPDIR="$WORK_DIR/ChiralCraft.AppDir"
  cp -a "$STAGE" "$APPDIR"
  cp -f "$ROOT_DIR/assets/chiralcraft.desktop" "$APPDIR/chiralcraft.desktop"
  cp -f "$ROOT_DIR/assets/chiralcraft-icon-128.png" "$APPDIR/chiralcraft.png"
  ln -sf chiralcraft.png "$APPDIR/.DirIcon"
  cat > "$APPDIR/AppRun" <<'EOF'
#!/bin/sh
HERE="$(dirname "$(readlink -f "$0")")"
exec "$HERE/usr/bin/chiralcraft" "$@"
EOF
  chmod 755 "$APPDIR/AppRun"
  ARCH="$APPIMAGE_ARCH" APPIMAGE_EXTRACT_AND_RUN=1 \
    "$APPIMAGETOOL" --no-appstream "$APPDIR" "$OUT_DIR/${PKG_BASE}.AppImage"
fi

echo "==> Linux packaging complete for ${ARCH}:"
ls -lh "$OUT_DIR/${PKG_BASE}"*
