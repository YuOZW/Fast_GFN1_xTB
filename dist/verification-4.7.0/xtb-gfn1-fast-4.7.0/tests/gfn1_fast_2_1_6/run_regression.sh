#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${TMPDIR:-/tmp}/gfn1_fast_216_policy_test"
rm -rf "$BUILD"
mkdir -p "$BUILD"
python3 "$ROOT/tests/gfn1_fast_2_1_6/test_source_guards.py"
gfortran -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan \
  -J"$BUILD" \
  "$ROOT/src/gfn1_fast_policy.f90" \
  "$ROOT/tests/gfn1_fast_2_1_6/test_spectral_policy.f90" \
  -o "$BUILD/test_spectral_policy"
"$BUILD/test_spectral_policy"
