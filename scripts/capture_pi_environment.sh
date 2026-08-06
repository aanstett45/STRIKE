#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
out="${1:-results/raspberry_pi/environment.txt}"
mkdir -p "$(dirname "$out")"

if command -v sha256sum >/dev/null 2>&1; then
  source_hash="$(sha256sum src/finger_spgemm.cpp | awk '{print $1}')"
else
  source_hash="$(shasum -a 256 src/finger_spgemm.cpp | awk '{print $1}')"
fi
commit="not available"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  commit="$(git rev-parse HEAD)"
fi
expected_governor="${PI_GOVERNOR:-performance}"

{
  echo "Measurement environment — Raspberry Pi"
  echo "======================================"
  echo
  echo "Capture time (UTC)"
  echo "------------------"
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  echo
  echo "Source identity"
  echo "---------------"
  echo "Git commit: $commit"
  echo "SHA-256 src/finger_spgemm.cpp: $source_hash"
  echo
  echo "Device model"
  echo "------------"
  tr -d '\0' < /proc/device-tree/model || true
  echo
  echo
  echo "Operating system"
  echo "----------------"
  cat /etc/os-release || true
  echo
  echo "Kernel (hostname omitted)"
  echo "-------------------------"
  uname -smr
  echo "Userland bits: $(getconf LONG_BIT 2>/dev/null || echo unknown)"
  echo
  echo "lscpu"
  echo "-----"
  lscpu
  echo
  echo "Memory"
  echo "------"
  free -h
  echo
  echo "g++ --version"
  echo "-------------"
  g++ --version
  echo
  echo "CPU governors"
  echo "-------------"
  echo "Expected CPU governor: $expected_governor"
  for f in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
    [ -r "$f" ] && printf '%s: %s\n' "$f" "$(cat "$f")"
  done
  echo
  echo "Throttling and temperature"
  echo "--------------------------"
  if command -v vcgencmd >/dev/null 2>&1; then
    vcgencmd get_throttled
    vcgencmd measure_temp
  else
    echo "vcgencmd unavailable"
  fi
  echo
  echo "Compilation command"
  echo "-------------------"
  echo 'g++ -O3 -march=native -fopenmp -std=c++17 src/finger_spgemm.cpp -o build/finger_spgemm'
  echo
  echo "Measurement conditions"
  echo "----------------------"
  echo "OMP_SCHEDULE=${OMP_SCHEDULE:-static}"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND:-close}"
  echo "OMP_PLACES=${OMP_PLACES:-cores}"
  echo "Power supply: ${PI_POWER_NOTES:-not provided}"
  echo "Cooling: ${PI_COOLING_NOTES:-not provided}"
  echo "Notes: ${MEASUREMENT_NOTES:-not provided}"
} > "$out"

printf 'Wrote %s\n' "$out"
