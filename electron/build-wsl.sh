#!/usr/bin/env bash
# Build SPLACE Desktop ZIPs from WSL.
# Must be run from the electron/ directory (or pass --win / --mac / --linux).
# Builds on the native Linux FS (/tmp) to avoid NTFS permission issues,
# then copies the ZIP back to electron/dist/.

set -e

PLATFORM="${1:---win}"   # default: --win

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="/tmp/splace-wsl-build"

echo "→ Copying project to Linux FS: $BUILD_DIR"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

# Copy electron/ and docs/ (needed by package.json "files" glob)
cp -r "$REPO_ROOT/electron" "$BUILD_DIR/electron"
cp -r "$REPO_ROOT/docs"     "$BUILD_DIR/docs"

echo "→ Installing npm dependencies"
cd "$BUILD_DIR/electron"
npm install --prefer-offline 2>/dev/null || npm install

echo "→ Building $PLATFORM ZIP"
CSC_IDENTITY_AUTO_DISCOVERY=false npx electron-builder "$PLATFORM"

echo "→ Copying output back to $SCRIPT_DIR/dist/"
mkdir -p "$SCRIPT_DIR/dist"
# Copy ZIPs, AppImages, and Windows portable EXEs
for ext in zip AppImage exe; do
    for f in "$BUILD_DIR/electron/dist/"*.$ext; do
        [ -f "$f" ] && cp "$f" "$SCRIPT_DIR/dist/"
    done
done

# Windows: also zip the portable exe for easy distribution
if [[ "$PLATFORM" == "--win" ]]; then
    cd "$SCRIPT_DIR/dist"
    for exe in *.exe; do
        [ -f "$exe" ] || continue
        zip "${exe%.exe}-windows-x64.zip" "$exe"
        echo "  Zipped: ${exe%.exe}-windows-x64.zip"
    done
    cd - > /dev/null
fi

echo ""
echo "Done! Output in $SCRIPT_DIR/dist/:"
ls -lh "$SCRIPT_DIR/dist/" | grep -E '\.(zip|AppImage|exe)' || ls -lh "$SCRIPT_DIR/dist/"
