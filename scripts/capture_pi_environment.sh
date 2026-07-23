#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
mkdir -p results/raspberry_pi
out=results/raspberry_pi/environment.txt

{
  echo "Measurement environment — Raspberry Pi"
  echo "======================================"
  echo
  echo "Known configuration"
  echo "-------------------"
  echo "Model: Raspberry Pi 4 Model B"
  echo "Memory: 4 GB"
  echo
  echo "Capture time"
  echo "------------"
  date
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
  echo "uname -a"
  echo "--------"
  uname -a
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
  for f in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
    [ -r "$f" ] && printf '%s: %s\n' "$f" "$(cat "$f")"
  done
  echo
  echo "Throttling status"
  echo "-----------------"
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
  echo "CPU governor required: performance"
  echo "OMP_SCHEDULE=static"
  echo "OMP_PROC_BIND=close"
  echo "OMP_PLACES=cores"
  echo "Power supply and cooling notes: EDIT THIS LINE BEFORE RELEASE"
} > "$out"

printf 'Wrote %s\n' "$out"
