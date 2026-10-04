#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEST="$ROOT/subprojects/mctc-lib"
MCTCLIB_COMMIT="${MCTCLIB_COMMIT:-77f65c6f2cf6330d05d0757ca173da097096780e}"
MCTCLIB_URL="${MCTCLIB_URL:-https://github.com/grimme-lab/mctc-lib.git}"

# xtb 6.7.1 was released on 2024-07-23 while its wrap requested mctc-lib
# revision=head. 77f65c6... was the mctc-lib main-branch HEAD immediately
# before that release date, so pinning this exact commit reproduces the
# dependency state much better than fetching today's HEAD.

if [[ -f "$DEST/CMakeLists.txt" ]]; then
  if [[ -d "$DEST/.git" ]]; then
    actual="$(git -C "$DEST" rev-parse HEAD)"
    if [[ "$actual" != "$MCTCLIB_COMMIT" ]]; then
      echo "ERROR: existing mctc-lib is not the pinned commit." >&2
      echo "Actual:   $actual" >&2
      echo "Expected: $MCTCLIB_COMMIT" >&2
      exit 3
    fi
  fi
  echo "Pinned/local mctc-lib source already present: $DEST"
  exit 0
fi

if [[ -e "$DEST" ]]; then
  echo "ERROR: $DEST exists but is not a usable mctc-lib source tree." >&2
  exit 2
fi

command -v git >/dev/null 2>&1 || {
  echo "ERROR: git is required to prepare the pinned dependency." >&2
  exit 2
}

mkdir -p "$DEST"
git -C "$DEST" init -q
git -C "$DEST" remote add origin "$MCTCLIB_URL"
echo "Fetching pinned mctc-lib commit ${MCTCLIB_COMMIT} ..."
if ! git -C "$DEST" fetch --depth 1 origin "$MCTCLIB_COMMIT"; then
  rm -rf "$DEST"
  exit 128
fi
git -C "$DEST" checkout -q --detach FETCH_HEAD
actual="$(git -C "$DEST" rev-parse HEAD)"
if [[ "$actual" != "$MCTCLIB_COMMIT" ]]; then
  echo "ERROR: fetched unexpected mctc-lib commit: $actual" >&2
  rm -rf "$DEST"
  exit 3
fi

echo "Pinned dependency ready: $DEST"
echo "Commit: $actual"
