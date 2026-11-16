#!/usr/bin/env bash
# Builds SPLACE for Linux (x64) and packs it as one portable file:
#   dist/SPLACE-<version>-x86_64.AppImage
# Runs on any recent distribution without installation (chmod +x, then run
# it). It carries the app, MAFFT (built from the official source), trimAl,
# IQ-TREE, RAxML-NG, ClipKIT and PhyKIT; the system provides GTK 3.
#
#   bash tool/package/appimage.sh
#
# Needs the Flutter SDK for Linux and: clang cmake ninja-build pkg-config
# libgtk-3-dev python3-venv gcc make curl. SKIP_BUILD=1 packs the last
# release build again; NO_PYTHON_TOOLS=1 leaves ClipKIT and PhyKIT out (the
# app then installs them with pip). MAFFT_VERSION picks the MAFFT release.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
version=$(sed -n "s/^const appVersion = '\(.*\)';/\1/p" lib/src/ui/common.dart)
mafft_version="${MAFFT_VERSION:-7.526}"
appimagetool_version="1.9.0"

[ "${SKIP_BUILD:-0}" = 1 ] || flutter build linux --release
bundle="${BUNDLE:-build/linux/x64/release/bundle}"
[ -x "$bundle/splace" ] || { echo "No release build in $bundle" >&2; exit 1; }
if [ "${NO_PYTHON_TOOLS:-0}" != 1 ] && [ ! -d "$bundle/tools/clipkit" ]; then
  bash "$root/tool/package/python_tools.sh" "$bundle"
fi

work="$root/build/appimage"
appdir="$work/SPLACE.AppDir"
rm -rf "$appdir"
mkdir -p "$appdir" "$work/cache"
cp -a "$bundle/." "$appdir/"

# MAFFT: the .deb takes it from the distribution, the AppImage carries its
# own (tools/mafft-linux, called through the tools/mafft launcher).
if [ ! -e "$appdir/tools/mafft" ]; then
  src="$work/cache/mafft-v$mafft_version"
  if [ ! -d "$src" ]; then
    curl -fsSL "https://gitlab.com/sysimm/mafft/-/archive/v$mafft_version/mafft-v$mafft_version.tar.gz" -o "$work/cache/mafft.tar.gz"
    tar -xzf "$work/cache/mafft.tar.gz" -C "$work/cache"
  fi
  make -C "$src/core" clean >/dev/null
  make -C "$src/core" PREFIX="$work/mafft-linux" -j"$(nproc)"
  make -C "$src/core" PREFIX="$work/mafft-linux" install
  mkdir -p "$appdir/tools"
  rm -rf "$appdir/tools/mafft-linux"
  cp -a "$work/mafft-linux" "$appdir/tools/mafft-linux"
  rm -rf "$appdir/tools/mafft-linux/share"
  cat > "$appdir/tools/mafft" <<'LAUNCHER'
#!/bin/sh
# Launcher for the bundled MAFFT: MAFFT needs MAFFT_BINARIES to find its
# helper programs.
here=$(cd "$(dirname "$0")" && pwd)
MAFFT_BINARIES="$here/mafft-linux/libexec/mafft"
export MAFFT_BINARIES
exec "$here/mafft-linux/bin/mafft" "$@"
LAUNCHER
fi

if [ -d "$appdir/tools" ]; then
  # Launcher scripts checked out on Windows can have CRLF line ends.
  for f in "$appdir/tools/"*; do
    [ -f "$f" ] && [ "$(head -c2 "$f")" = '#!' ] && sed -i 's/\r$//' "$f"
  done
  find "$appdir/tools" -type f ! -name '*.txt' -exec chmod 755 {} +
fi
chmod 755 "$appdir/splace"

cat > "$appdir/AppRun" <<'APPRUN'
#!/bin/sh
here=$(dirname "$(readlink -f "$0")")
exec "$here/splace" "$@"
APPRUN
chmod 755 "$appdir/AppRun"
cp assets/images/app_icon.png "$appdir/splace.png"
ln -sf splace.png "$appdir/.DirIcon"
cat > "$appdir/splace.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=SPLACE
Comment=SPLit, Align and ConcatenatE: phylogenomic data preparation
Exec=splace %F
Icon=splace
Categories=Science;Biology;
Terminal=false
StartupWMClass=splace
X-AppImage-Version=$version
DESKTOP

tool="$work/cache/appimagetool-$appimagetool_version-x86_64.AppImage"
if [ ! -x "$tool" ]; then
  curl -fsSL "https://github.com/AppImage/appimagetool/releases/download/$appimagetool_version/appimagetool-x86_64.AppImage" -o "$tool"
  chmod 755 "$tool"
fi
mkdir -p dist
out="dist/SPLACE-$version-x86_64.AppImage"
# Extract-and-run: works where FUSE is not available (containers, CI).
APPIMAGE_EXTRACT_AND_RUN=1 ARCH=x86_64 "$tool" --no-appstream "$appdir" "$out"
echo "AppImage: $root/$out"
