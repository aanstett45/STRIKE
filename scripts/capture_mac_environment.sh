#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p results/macbook
out=results/macbook/environment.txt

{
  echo "Measurement environment — MacBook"
  echo "================================="
  echo
  echo "Known configuration"
  echo "-------------------"
  echo "Model: MacBook Pro"
  echo "SoC: Apple M4 Pro"
  echo "Memory: 24 GB"
  echo "Operating system: macOS 15.7.4"
  echo "Build: 24G517"
  echo
  echo "Capture time"
  echo "------------"
  date
  echo
  echo "system_profiler SPHardwareDataType SPSoftwareDataType"
  echo "----------------------------------------------------"
  system_profiler SPHardwareDataType SPSoftwareDataType
  echo
  echo "uname -a"
  echo "--------"
  uname -a
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
  echo "Power: connected to AC power"
  echo "Sleep: disabled during measurements"
  echo "Default schedule: OMP_SCHEDULE=static"
  echo "Additional notes: EDIT THIS LINE BEFORE RELEASE"
} > "$out"

printf 'Wrote %s\n' "$out"
