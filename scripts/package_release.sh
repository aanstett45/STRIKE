#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

bash scripts/prepare_release.sh
out="${1:-../STRIKE-v1.0.1.zip}"
rm -f "$out"
zip -qr "$out" . \
  -x '.git/*' 'build/*' '.benchmark_backups/*' '*/__pycache__/*' '*.pyc' '*.DS_Store' '*.zip'
echo "Wrote $out"
