#!/usr/bin/env bash
# Builds trimAl for Linux from the official source at a given tag (default
# v1.5.0, the same revision as the Windows trimal.exe shipped with SPLACE),
# statically linked so it runs on any x86-64 Linux, tests the case that
# crashes trimAl v1.5.1 (-backtrans with -colnumbering) and, when it
# passes, replaces GUI/tools/linux/trimal.
#
#   bash tool/package/build_trimal_linux.sh [tag]
#
# Needs git, g++ and make (sudo apt install git build-essential).
set -euo pipefail
tag="${1:-v1.5.0}"
root="$(cd "$(dirname "$0")/../.." && pwd)"
bin="$root/tools/linux"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

git clone -q --depth 1 --branch "$tag" https://github.com/inab/trimal.git "$work/src"
make -s -C "$work/src/source" FLAGS="-O2 -static -w" trimal
t="$work/src/source/trimal"
"$t" --version

# A codon alignment (2 codons broken by the gap of c) and its translation.
cat > "$work/prot.fa" <<'FA'
>a
MAYSLGM
>b
MAYSLGM
>c
MAN-LGM
FA
cat > "$work/cds.fa" <<'FA'
>a
ATGGCTTATTCACTGGGAATG
>b
ATGGCCTATTCACTGGGAATG
>c
ATGGCTAATCTGGGAATG
FA
"$t" -in "$work/prot.fa" -out "$work/out.fa" -backtrans "$work/cds.fa" -gappyout -ignorestopcodon -colnumbering > "$work/map.txt"
grep -q '#ColumnsMap' "$work/map.txt" || { echo "trimAl $tag printed no column map" >&2; exit 1; }
len=$(awk '!/^>/{printf "%s", $0; exit}' "$work/out.fa" | wc -c)
[ $((len % 3)) -eq 0 ] || { echo "trimAl $tag: back-translated length $len is not a multiple of 3" >&2; exit 1; }
echo "trimAl $tag: -backtrans with -colnumbering OK ($len bp of whole codons)"

install -m 755 "$work/src/source/trimal" "$bin/trimal"
echo "Installed in $bin"
