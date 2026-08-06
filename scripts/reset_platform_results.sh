#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 || ( "$1" != "macbook" && "$1" != "raspberry_pi" ) ]]; then
  echo "Usage: $0 macbook|raspberry_pi" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
platform="$1"
dir="results/$platform"

if [[ -d "$dir" && -n "$(find "$dir" -type f ! -name '.gitkeep' -print -quit)" ]]; then
  stamp="$(date -u +%Y%m%dT%H%M%SZ)"
  backup=".benchmark_backups/${platform}-${stamp}"
  mkdir -p .benchmark_backups
  mv "$dir" "$backup"
  echo "Previous $platform results moved to $backup"
fi

mkdir -p "$dir/threads" "$dir/scheduling"
touch "$dir/.gitkeep" "$dir/threads/.gitkeep" "$dir/scheduling/.gitkeep"
