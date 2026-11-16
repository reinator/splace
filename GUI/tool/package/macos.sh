#!/usr/bin/env bash
# Builds SPLACE for macOS (universal: Apple silicon and Intel) and packs it
# as dist/SPLACE-<version>-macos-universal.dmg, with MAFFT, trimAl,
# IQ-TREE and RAxML-NG inside the app. Run on a Mac with Xcode:
#
#   bash tool/package/macos.sh
#
# SKIP_BUILD=1 packs the last release build again. Without an Apple
# Developer ID the app is signed ad hoc: on another Mac, open it the first
# time with right-click > Open (or: xattr -dr com.apple.quarantine /Applications/SPLACE.app).
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
version=$(sed -n "s/^const appVersion = '\(.*\)';/\1/p" lib/src/ui/common.dart)

[ "${SKIP_BUILD:-0}" = 1 ] || flutter build macos --release
app="build/macos/Build/Products/Release/splace.app"
[ -d "$app" ] || { echo "No release build in $app" >&2; exit 1; }
if /usr/libexec/PlistBuddy -c 'Print :com.apple.security.app-sandbox' macos/Runner/Release.entitlements 2>/dev/null | grep -q true; then
  echo "WARNING: Release.entitlements has app-sandbox = true: the bundled tools and NCBI downloads will not run." >&2
fi

stage="$(mktemp -d)"
trap 'rm -rf "$stage"' EXIT
cp -R "$app" "$stage/SPLACE.app"
ln -s /Applications "$stage/Applications"
# Apple silicon refuses unsigned code: sign the bundled tools, then the app
# (ad hoc, keeping its entitlements).
find "$stage/SPLACE.app/Contents/Resources/tools" -type f -perm -u+x -print0 2>/dev/null |
  while IFS= read -r -d '' f; do codesign --force --sign - "$f" 2>/dev/null || true; done
codesign --force --sign - --preserve-metadata=entitlements "$stage/SPLACE.app"
xattr -dr com.apple.quarantine "$stage/SPLACE.app" 2>/dev/null || true

mkdir -p dist
hdiutil create -volname "SPLACE $version" -srcfolder "$stage" -ov -format UDZO "dist/SPLACE-$version-macos-universal.dmg"
echo "Disk image: $root/dist/SPLACE-$version-macos-universal.dmg"
