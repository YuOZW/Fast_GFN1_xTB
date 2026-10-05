#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${TMPDIR:-/tmp}/gfn1_fast_226_policy_test"
rm -rf "$BUILD"
mkdir -p "$BUILD"
python3 "$ROOT/tests/gfn1_fast_2_2_6/test_source_guards.py"
python3 "$ROOT/tests/gfn1_fast_2_2_6/test_compact_density_math.py"
gfortran -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan \
  -J"$BUILD" \
  "$ROOT/src/gfn1_fast_policy.f90" \
  "$ROOT/tests/gfn1_fast_2_2_6/test_partial_dispatch_policy.f90" \
  -o "$BUILD/test_partial_dispatch_policy"
"$BUILD/test_partial_dispatch_policy"
