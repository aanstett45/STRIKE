#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
out="${1:-results/macbook/environment.txt}"
mkdir -p "$(dirname "$out")"

source_hash="$(shasum -a 256 src/finger_spgemm.cpp | awk '{print $1}')"
commit="not available"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  commit="$(git rev-parse HEAD)"
fi

{
  echo "Measurement environment — MacBook"
  echo "================================="
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
  echo "Hardware and operating system (privacy-filtered)"
  echo "------------------------------------------------"
  system_profiler SPHardwareDataType SPSoftwareDataType | sed -E \
    -e '/Serial Number:/d' \
    -e '/Hardware UUID:/d' \
    -e '/Provisioning UDID:/d' \
    -e '/Computer Name:/d' \
    -e '/Local Host Name:/d' \
    -e '/Host Name:/d' \
    -e '/User Name:/d' \
    -e '/Activation Lock Status:/d'
  echo
  echo "Kernel (hostname omitted)"
  echo "-------------------------"
  uname -smr
  echo
  echo "clang++ --version"
  echo "-----------------"
  clang++ --version
  echo
  echo "brew list --versions libomp"
  echo "---------------------------"
  brew list --versions libomp
  echo
  echo "Compilation command"
  echo "-------------------"
  echo 'clang++ -O3 -march=native -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp -std=c++17 src/finger_spgemm.cpp -o build/finger_spgemm'
  echo
  echo "Measurement conditions"
  echo "----------------------"
  echo "Power: ${MAC_POWER_STATE:-connected to AC power}"
  echo "Sleep: ${MAC_SLEEP_STATE:-disabled during measurements}"
  echo "OMP_SCHEDULE=${OMP_SCHEDULE:-static}"
  echo "Notes: ${MEASUREMENT_NOTES:-not provided}"
} > "$out"

printf 'Wrote %s\n' "$out"
