#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${1:-$ROOT/build-gfn1-fast-2.2.6-debug}"

# This debug route intentionally validates the native GFN1-xTB path only.
# tblite and CPCM-X are independent optional backends and are disabled so an
# old xtb-6.7.1 tree never downloads today's incompatible dependency HEADs.
# Unit tests are disabled here because test-drive is not needed to compile or
# execute the GFN1 Energy/Gradient path; the 2.2.6 standalone regressions are
# run separately by run_regression.sh.
#
# Dependency lookup deliberately omits "fetch".  Supply mctc-lib either as an
# installed package or as ROOT/subprojects/mctc-lib (prepare_debug_deps.sh
# creates the latter at a pinned tag when networking is available).

if [[ ! -f "$ROOT/subprojects/mctc-lib/CMakeLists.txt" ]]; then
  echo "NOTE: no vendored mctc-lib found at $ROOT/subprojects/mctc-lib." >&2
  echo "      CMake will try installed mctc-lib locations only; it will NOT fetch HEAD." >&2
fi

CMAKE_EXTRA=()
if [[ -n "${MCTCLIB_DIR:-}" ]]; then
  CMAKE_EXTRA+=("-Dmctc-lib_DIR=${MCTCLIB_DIR}")
fi

rm -rf "$BUILD"
cmake -S "$ROOT" -B "$BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DWITH_TBLITE=OFF \
  -DWITH_CPCMX=OFF \
  -DWITH_TESTS=OFF \
  -DWITH_JSON=OFF \
  '-DMCTCLIB_FIND_METHOD=subproject' \
  '-DCMAKE_Fortran_FLAGS_DEBUG=-O0 -g3 -fcheck=all -finit-real=snan -finit-integer=-999999 -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-omit-frame-pointer' \
  "${CMAKE_EXTRA[@]}"

cmake --build "$BUILD" -j"${CMAKE_BUILD_PARALLEL_LEVEL:-2}"

cat <<'MSG'
GFN1-fast 2.2.6 native-GFN1 full Debug build completed with:
  - bounds/runtime checks
  - signaling-NaN REAL initialization
  - sentinel INTEGER initialization
  - floating-point traps
  - backtraces/frame pointers
  - no unpinned dependency fetches
Run tests/gfn1_fast_2_2_6/run_regression.sh and then the normal Energy/Gradient
molecule regression set with the Debug fast-gfn1-xtb executable.
MSG
