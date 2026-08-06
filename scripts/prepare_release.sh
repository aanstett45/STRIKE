#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ -d .benchmark_backups && -n "$(find .benchmark_backups -type f -print -quit)" ]]; then
  echo "Remove or move .benchmark_backups before packaging a release." >&2
  exit 1
fi

python3 scripts/validate_results.py --platform all --compare-platforms --pi-governor "${PI_GOVERNOR:-performance}"
python3 scripts/make_figure1.py results/macbook/fullspgemm.csv > generated/figure1_addplots.tex

if grep -RIn --exclude='CHANGELOG.md' --exclude='MANIFEST.sha256' --exclude='prepare_release.sh' \
  --exclude-dir='.git' --exclude-dir='build' --exclude-dir='.benchmark_backups' \
  'EDIT THIS LINE BEFORE RELEASE' .; then
  echo "Release placeholder remains." >&2
  exit 1
fi

bash scripts/regenerate_manifest.sh

echo "Release checks passed. Review environment notes and article values before tagging v1.0.1."
