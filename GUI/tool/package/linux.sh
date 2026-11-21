#!/usr/bin/env bash
# Builds SPLACE for Linux (x64) and packs it as dist/splace_<version>_amd64.deb
# (app in /opt/splace with IQ-TREE, RAxML-NG, trimAl, ClipKIT and PhyKIT;
# MAFFT comes from the distribution's package, a dependency of the .deb).
#
#   bash tool/package/linux.sh
#
# Needs the Flutter SDK for Linux and: clang cmake ninja-build pkg-config
# libgtk-3-dev dpkg-dev python3-venv. SKIP_BUILD=1 packs the last release
# build again; NO_PYTHON_TOOLS=1 leaves ClipKIT and PhyKIT out (the app
# then installs them with pip).
# Install the result with:  sudo apt install ./dist/splace_<version>_amd64.deb
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
version=$(sed -n "s/^const appVersion = '\(.*\)';/\1/p" lib/src/ui/common.dart)

[ "${SKIP_BUILD:-0}" = 1 ] || flutter build linux --release
bundle="${BUNDLE:-build/linux/x64/release/bundle}"
[ -x "$bundle/splace" ] || { echo "No release build in $bundle" >&2; exit 1; }

# ClipKIT and PhyKIT as standalone programs (PyInstaller, one folder each:
# tools/<name>/<name>), so the users need no Python.
if [ "${NO_PYTHON_TOOLS:-0}" != 1 ]; then
  bash "$root/tool/package/python_tools.sh" "$bundle"
fi

stage="$(mktemp -d)"
trap 'rm -rf "$stage"' EXIT
pkg="$stage/splace_${version}_amd64"
mkdir -p "$pkg/DEBIAN" "$pkg/opt/splace" "$pkg/usr/bin" "$pkg/usr/share/applications" "$pkg/usr/share/icons/hicolor/512x512/apps" "$pkg/usr/share/mime/packages"
cp -a "$bundle/." "$pkg/opt/splace/"
if [ -d "$pkg/opt/splace/tools" ]; then
  # Launcher scripts checked out on Windows can have CRLF line ends.
  for f in "$pkg/opt/splace/tools/"*; do
    [ -f "$f" ] && [ "$(head -c2 "$f")" = '#!' ] && sed -i 's/\r$//' "$f"
  done
  find "$pkg/opt/splace/tools" -type f ! -name '*.txt' -exec chmod 755 {} +
fi
chmod 755 "$pkg/opt/splace/splace"
ln -s /opt/splace/splace "$pkg/usr/bin/splace"
cp assets/images/app_icon.png "$pkg/usr/share/icons/hicolor/512x512/apps/splace.png"
cat > "$pkg/usr/share/applications/splace.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=SPLACE
Comment=SPLit, Align and ConcatenatE: phylogenomic data preparation
Exec=/opt/splace/splace %F
Icon=splace
Categories=Science;Biology;
MimeType=application/x-splace-project;
Terminal=false
StartupWMClass=splace
DESKTOP
# Project files (.splc) open in SPLACE (dpkg triggers update the MIME and
# desktop databases).
cat > "$pkg/usr/share/mime/packages/splace.xml" <<MIME
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-splace-project">
    <comment>SPLACE project</comment>
    <sub-class-of type="application/zip"/>
    <glob pattern="*.splc" weight="60"/>
  </mime-type>
</mime-info>
MIME
cat > "$pkg/DEBIAN/control" <<CONTROL
Package: splace
Version: $version
Section: science
Priority: optional
Architecture: amd64
Maintainer: Luan Rabelo <lprabelo.itv@gmail.com>
Installed-Size: $(du -sk "$pkg" | cut -f1)
Depends: libgtk-3-0, libglib2.0-0, mafft
Recommends: python3, python3-pip
Homepage: https://github.com/luanrabelo/splace
Description: SPLit, Align and ConcatenatE
 Phylogenomic data preparation: GenBank and FASTA input, marker
 extraction, MAFFT alignment, trimming, saturation test, supermatrix
 and trees with IQ-TREE and RAxML-NG (bundled).
CONTROL
mkdir -p dist
dpkg-deb --build --root-owner-group "$pkg" "dist/splace_${version}_amd64.deb"
echo "Package: $root/dist/splace_${version}_amd64.deb"
