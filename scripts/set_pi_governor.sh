#!/usr/bin/env bash
set -euo pipefail

target="${1:-performance}"
case "$target" in
  performance|ondemand|powersave|conservative|schedutil) ;;
  *) echo "Unsupported governor name: $target" >&2; exit 2 ;;
esac

mapfile -t files < <(compgen -G '/sys/devices/system/cpu/cpu*/cpufreq/scaling_governor' || true)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No CPU governor controls found." >&2
  exit 1
fi

write_one() {
  local f="$1"
  if [[ -w "$f" ]]; then
    printf '%s\n' "$target" > "$f"
  else
    printf '%s\n' "$target" | sudo tee "$f" >/dev/null
  fi
}

for f in "${files[@]}"; do write_one "$f"; done

for f in "${files[@]}"; do
  actual="$(cat "$f")"
  if [[ "$actual" != "$target" ]]; then
    echo "Governor verification failed for $f: expected $target, got $actual" >&2
    exit 1
  fi
done

echo "All CPU governors set to $target."
