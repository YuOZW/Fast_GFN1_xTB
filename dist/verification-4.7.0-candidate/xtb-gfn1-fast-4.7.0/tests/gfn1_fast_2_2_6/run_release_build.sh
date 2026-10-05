#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${1:-$ROOT/build-gfn1-fast-2.2.6-release}"

if [[ ! -f "$ROOT/subprojects/mctc-lib/CMakeLists.txt" ]]; then
  echo "ERROR: pinned/local mctc-lib source is required at $ROOT/subprojects/mctc-lib" >&2
  echo "Run tests/gfn1_fast_2_2_6/prepare_debug_deps.sh first." >&2
  exit 2
fi

rm -rf "$BUILD"
cmake -S "$ROOT" -B "$BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_Fortran_FLAGS_RELEASE=-O3 -march=native -mtune=native' \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF \
  -DMCTCLIB_FIND_METHOD=subproject \
  -DWITH_TBLITE=OFF \
  -DWITH_CPCMX=OFF \
  -DWITH_TESTS=OFF \
  -DWITH_JSON=OFF
cmake --build "$BUILD" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}"

echo "GFN1-fast 2.2.6 Release executable: $BUILD/xtb"
