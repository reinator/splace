#!/bin/sh
# Copies the external tools (MAFFT, trimAl, IQ-TREE, RAxML-NG) into
# SPLACE.app/Contents/Resources/tools, where the app looks for them first.
# SPLACE_TOOLS_DIR overrides the source folder (default: the repository's
# electron/bin/mac). Run by Xcode as a build phase of the Runner target.
set -e
src="${SPLACE_TOOLS_DIR:-$SRCROOT/../../electron/bin/mac}"
[ -d "$src" ] || { echo "SPLACE: no tools in $src (skipped)"; exit 0; }
dest="$TARGET_BUILD_DIR/$UNLOCALIZED_RESOURCES_FOLDER_PATH/tools"
mkdir -p "$dest"
rsync -a --delete --exclude seqkit --exclude readal --exclude statal "$src/" "$dest/"
chmod +x "$dest"/* "$dest"/mafft-mac/mafftdir/bin/* "$dest"/mafft-mac/mafftdir/libexec/* 2>/dev/null || true
# Downloaded binaries keep the quarantine flag; the app's copies must not.
xattr -dr com.apple.quarantine "$dest" 2>/dev/null || true
echo "SPLACE: tools copied to $dest"
