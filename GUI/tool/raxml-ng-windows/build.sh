#!/bin/sh
# Builds raxml-ng.exe (RAxML-NG 2.0.4, 64-bit Windows) on Linux/WSL with the
# llvm-mingw cross compiler, from the official source and the patches here:
#   patch_source.py  Windows fixes in RAxML-NG/coraxlib (system calls, a
#                    random index computed with RAND_MAX, static coraxlib)
#   win_compat.h     asprintf, glibc-compatible rand()/srand() (same seed,
#                    same random numbers as on Linux), aligned allocation
#   splace_alloc.c   one allocator for the program: coraxlib frees some
#                    aligned blocks with free(), which the Windows CRT rejects
# Needs: curl, unzip, python3, cmake and ninja (e.g. pip install cmake ninja).
# Output: ./out/raxml-ng.exe (depends only on Windows system DLLs).
set -e
here=$(cd "$(dirname "$0")" && pwd)
work=${WORK:-$here/work}
mkdir -p "$work" && cd "$work"
LLVM_MINGW=${LLVM_MINGW:-llvm-mingw-20261104-ucrt-ubuntu-22.04-x86_64}
if [ ! -d "$LLVM_MINGW" ]; then
  curl -sL -o tc.tar.xz "https://github.com/mstorsjo/llvm-mingw/releases/download/${LLVM_MINGW#llvm-mingw-}/$LLVM_MINGW.tar.xz" 2>/dev/null ||
    curl -sL -o tc.tar.xz "https://github.com/mstorsjo/llvm-mingw/releases/download/20261104/$LLVM_MINGW.tar.xz"
  tar -xJf tc.tar.xz
fi
TC="$work/$LLVM_MINGW"
rm -rf src && mkdir src && cd src
curl -sL -o src.zip https://github.com/amkozlov/raxml-ng/releases/download/2.0.4/raxml-ng_v2.0.4_source.zip
unzip -q src.zip && rm src.zip
python3 "$here/patch_source.py" "$work/src"
cp "$here/win_compat.h" "$work/src/"
cd "$work"
"$TC/bin/x86_64-w64-mingw32-clang" -O2 -c "$here/splace_alloc.c" -o splace_alloc.o
cat > mingw.cmake <<CM
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)
set(CMAKE_C_COMPILER $TC/bin/x86_64-w64-mingw32-clang)
set(CMAKE_CXX_COMPILER $TC/bin/x86_64-w64-mingw32-clang++)
set(CMAKE_RC_COMPILER $TC/bin/x86_64-w64-mingw32-windres)
set(CMAKE_FIND_ROOT_PATH $TC/x86_64-w64-mingw32)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
CM
rm -rf build && mkdir build && cd build
FLAGS="-D_USE_MATH_DEFINES -include $work/src/win_compat.h"
cmake -G Ninja -DCMAKE_TOOLCHAIN_FILE=../mingw.cmake -DSTATIC_BUILD=ON -DPORTABLE_BUILD=ON -DCMAKE_BUILD_TYPE=Release \
  "-DCMAKE_C_FLAGS=$FLAGS" "-DCMAKE_CXX_FLAGS=$FLAGS" "-DCMAKE_EXE_LINKER_FLAGS=-static $work/splace_alloc.o -lpsapi" ../src
ninja
mkdir -p "$here/out" && cp bin/raxml-ng.exe "$here/out/"
echo "Built $here/out/raxml-ng.exe"
