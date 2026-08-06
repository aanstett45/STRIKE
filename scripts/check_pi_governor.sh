#!/usr/bin/env bash
set -euo pipefail

target="${1:-performance}"
mapfile -t files < <(compgen -G '/sys/devices/system/cpu/cpu*/cpufreq/scaling_governor' || true)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No CPU governor controls found." >&2
  exit 1
fi

bad=0
for f in "${files[@]}"; do
  actual="$(cat "$f")"
  printf '%s: %s\n' "$f" "$actual"
  [[ "$actual" == "$target" ]] || bad=1
done

if [[ $bad -ne 0 ]]; then
  echo "Expected governor '$target'. Run: sudo bash scripts/set_pi_governor.sh '$target'" >&2
  exit 1
fi
