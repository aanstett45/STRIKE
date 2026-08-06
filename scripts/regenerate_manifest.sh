#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
tmp="MANIFEST.sha256.tmp"
: > "$tmp"

hash_file() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$f"
  else
    shasum -a 256 "$f"
  fi
}

while IFS= read -r -d '' file; do
  rel="${file#./}"
  digest="$(hash_file "$rel" | awk '{print $1}')"
  printf '%s  %s\n' "$digest" "$rel" >> "$tmp"
done < <(find . -type f \
  ! -path './.git/*' \
  ! -path './build/*' \
  ! -path './.benchmark_backups/*' \
  ! -path '*/__pycache__/*' \
  ! -name '.DS_Store' \
  ! -name 'MANIFEST.sha256' \
  ! -name 'MANIFEST.sha256.tmp' \
  ! -name '*.zip' \
  ! -name '*.pyc' \
  -print0 | sort -z)

mv "$tmp" MANIFEST.sha256
echo "Wrote MANIFEST.sha256"
