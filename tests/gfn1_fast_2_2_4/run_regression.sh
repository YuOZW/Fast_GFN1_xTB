#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${TMPDIR:-/tmp}/gfn1_fast_224_policy_test"
rm -rf "$BUILD"
mkdir -p "$BUILD"
python3 "$ROOT/tests/gfn1_fast_2_2_4/test_source_guards.py"
gfortran -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan \
  -J"$BUILD" \
  "$ROOT/src/gfn1_fast_policy.f90" \
  "$ROOT/tests/gfn1_fast_2_2_4/test_partial_dispatch_policy.f90" \
  -o "$BUILD/test_partial_dispatch_policy"
"$BUILD/test_partial_dispatch_policy"
