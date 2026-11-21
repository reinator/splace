#!/usr/bin/env bash
# Turns ClipKIT, PhyKIT and GMYC into standalone programs (PyInstaller, one folder
# each: <bundle>/tools/<name>/<name>), so the users need no Python. Used by
# linux.sh and appimage.sh:
#
#   bash tool/package/python_tools.sh build/linux/x64/release/bundle
#
# Needs python3 with venv. The virtual environment is kept in build/pytools.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
bundle="$(cd "${1:?usage: python_tools.sh <bundle folder>}" && pwd)"

venv="$root/build/pytools"
[ -x "$venv/bin/python" ] || python3 -m venv "$venv"
"$venv/bin/pip" install --upgrade pip pyinstaller clipkit phykit
mkdir -p "$bundle/tools"
for t in clipkit:clipkit.clipkit phykit:phykit.phykit; do
  name="${t%%:*}" module="${t#*:}"
  printf 'import sys\nfrom %s import main\nsys.exit(main())\n' "$module" > "$venv/${name}_entry.py"
  "$venv/bin/pyinstaller" --noconfirm --onedir --name "$name" --collect-all "$name" --collect-data Bio \
    --distpath "$bundle/tools" --workpath "$venv/work" --specpath "$venv" "$venv/${name}_entry.py"
done
# GMYC (pGMYC, the iTaxoTools Python 3 port in tool/gmyc, with SPLACE's
# command line): ete3 and scipy only, no Qt or matplotlib.
"$venv/bin/pip" install --upgrade ete3 scipy six
"$venv/bin/pyinstaller" --noconfirm --onedir --name gmyc --paths "$root/tool/gmyc" --hidden-import GMYC \
  --collect-submodules ete3 --exclude-module PyQt5 --exclude-module matplotlib \
  --distpath "$bundle/tools" --workpath "$venv/work" --specpath "$venv" "$root/tool/gmyc/splace_gmyc.py"
