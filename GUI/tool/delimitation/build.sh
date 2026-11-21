#!/bin/sh
# Builds the command-line programs of ASAP (Puillandre, Brouillet & Achaz
# 2021) and ABGD (Puillandre et al. 2012) for SPLACE's species delimitation,
# from the C sources kept by iTaxoTools (GPL-3.0), into GUI/tools/<platform>:
#   sh tool/delimitation/build.sh linux     static x86-64 binaries (gcc)
#   sh tool/delimitation/build.sh windows   asap.exe/abgd.exe with llvm-mingw
#   sh tool/delimitation/build.sh mac       universal binaries (clang, on a Mac)
# The ASAP sources carry iTaxoTools' Python output wrapper (wrapio.h); the
# command-line build replaces it with an empty header (plain stdio).
# Needs curl and tar (+ gcc, or clang on macOS; windows downloads llvm-mingw).
set -e
target=${1:-linux}
here=$(cd "$(dirname "$0")" && pwd)
gui=$(cd "$here/../.." && pwd)
work=${WORK:-$here/work}
ASAP_REV=7cd38ba422513ee4ab0f4b4672a661edee89dd45
ABGD_REV=7452ec25f5523479f0850dc97d32e59f3c1b2d3f
mkdir -p "$work" && cd "$work"
fetch() { # repo rev dir
  if [ ! -d "$3" ]; then
    curl -sSL -o "$3.tar.gz" "https://codeload.github.com/iTaxoTools/$1/tar.gz/$2"
    mkdir -p "$3" && tar -xzf "$3.tar.gz" -C "$3" --strip-components=1 && rm "$3.tar.gz"
  fi
}
fetch ASAPy "$ASAP_REV" asapy
fetch ABGDpy "$ABGD_REV" abgdpy
rm -rf asap abgd && cp -r asapy/src/asap asap && cp -r abgdpy/src/abgd abgd
echo '/* SPLACE: plain stdio (no Python wrapper) */' > asap/wrapio.h
# Old C: functions used before their declaration (an error in recent clang).
COMPAT="-Wno-error=implicit-function-declaration -Wno-error=int-conversion -Wno-error=incompatible-function-pointer-types"
ASAP_SRC="oldfns.c asap.c asap_common.c asap_core.c gdtosvg.c draw.c"
ABGD_SRC="abgdCore.c bionjcabgd.c main_abgd.c"

case "$target" in
  linux)
    out="$gui/tools/linux"; exe=""
    build() { gcc -O2 -w -static -DLINUX "$@" -lm; } ;;
  mac)
    out="$gui/tools/mac"; exe=""
    build() { clang -O2 -w $COMPAT -arch arm64 -arch x86_64 -mmacosx-version-min=11.0 -DMACOSX "$@" -lm; } ;;
  windows)
    out="$gui/tools/win"; exe=".exe"
    # llvm-mingw release (LLVM_MINGW_TAG=latest takes the newest one).
    tag=${LLVM_MINGW_TAG:-20261117}
    if [ "$tag" = latest ]; then
      tag=$(curl -sSL https://api.github.com/repos/mstorsjo/llvm-mingw/releases/latest | sed -n 's/.*"tag_name": *"\([^"]*\)".*/\1/p')
    fi
    LLVM_MINGW=llvm-mingw-$tag-ucrt-ubuntu-22.04-x86_64
    if [ ! -d "$LLVM_MINGW" ]; then
      curl -fsSL -o tc.tar.xz "https://github.com/mstorsjo/llvm-mingw/releases/download/$tag/$LLVM_MINGW.tar.xz"
      tar -xJf tc.tar.xz && rm tc.tar.xz
    fi
    CC="$work/$LLVM_MINGW/bin/x86_64-w64-mingw32-clang"
    build() { "$CC" -O2 -w $COMPAT -static -include "$here/splace_compat.h" "$@" "$here/splace_compat.c" -lm; } ;;
  *) echo "usage: $0 linux|windows|mac" >&2; exit 1 ;;
esac

mkdir -p "$out"
(cd asap && build -o "$out/asap$exe" $ASAP_SRC)
(cd abgd && build -o "$out/abgd$exe" $ABGD_SRC)
cp asapy/LICENSE.txt "$out/asap-abgd.LICENSE.txt"
echo "Built $out/asap$exe and $out/abgd$exe"
