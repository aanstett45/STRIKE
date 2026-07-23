#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

for f in data/roadNet-CA.txt data/email-Enron.txt; do
  if [[ ! -f "$f" ]]; then
    echo "Missing required dataset: $f" >&2
    echo "Run: bash scripts/download_snap_data.sh" >&2
    exit 1
  fi
done

mkdir -p build results/raspberry_pi/threads

g++ -O3 -march=native -fopenmp -std=c++17 \
  src/finger_spgemm.cpp -o build/finger_spgemm

bash scripts/capture_pi_environment.sh

export OMP_SCHEDULE=static
export OMP_PROC_BIND=close
export OMP_PLACES=cores

./build/finger_spgemm verify \
  2>&1 | tee results/raspberry_pi/verify.log

for k in 1 2 4; do
  ./build/finger_spgemm triangles \
    "results/raspberry_pi/threads/tri_T${k}.csv" \
    --runs 5 --threads "$k" \
    data/roadNet-CA.txt data/email-Enron.txt \
    gen:band:20000:16 gen:ba:20000:8 \
    2>&1 | tee "results/raspberry_pi/threads/tri_T${k}.log"
done

echo "Raspberry Pi campaign completed. Review logs and environment metadata before release."
