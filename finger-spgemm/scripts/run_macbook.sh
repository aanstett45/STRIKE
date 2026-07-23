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

mkdir -p build results/macbook/threads results/macbook/scheduling generated

clang++ -O3 -march=native \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -std=c++17 src/finger_spgemm.cpp \
  -o build/finger_spgemm

bash scripts/capture_mac_environment.sh

OMP_SCHEDULE=static ./build/finger_spgemm verify \
  2>&1 | tee results/macbook/verify.log

for k in 1 2 4 8 12; do
  OMP_SCHEDULE=static ./build/finger_spgemm triangles \
    "results/macbook/threads/tri_T${k}.csv" \
    --runs 5 --threads "$k" \
    data/roadNet-CA.txt data/email-Enron.txt \
    gen:band:20000:16 gen:ba:20000:8 \
    2>&1 | tee "results/macbook/threads/tri_T${k}.log"
done

OMP_SCHEDULE=static ./build/finger_spgemm triangles \
  results/macbook/scheduling/static.csv \
  --runs 5 --threads 8 \
  data/email-Enron.txt gen:ba:20000:8 \
  2>&1 | tee results/macbook/scheduling/static.log

OMP_SCHEDULE=dynamic,64 ./build/finger_spgemm triangles \
  results/macbook/scheduling/dynamic64.csv \
  --runs 5 --threads 8 \
  data/email-Enron.txt gen:ba:20000:8 \
  2>&1 | tee results/macbook/scheduling/dynamic64.log

OMP_SCHEDULE=guided ./build/finger_spgemm triangles \
  results/macbook/scheduling/guided.csv \
  --runs 5 --threads 8 \
  data/email-Enron.txt gen:ba:20000:8 \
  2>&1 | tee results/macbook/scheduling/guided.log

OMP_SCHEDULE=static ./build/finger_spgemm fullspgemm \
  results/macbook/fullspgemm.csv 2048 3 3 \
  2>&1 | tee results/macbook/fullspgemm.log

OMP_SCHEDULE=static ./build/finger_spgemm fullspgemm \
  results/macbook/full_wide.csv 4096 1 3 \
  2>&1 | tee results/macbook/full_wide.log

python3 scripts/make_figure1.py results/macbook/fullspgemm.csv \
  | tee generated/figure1_addplots.tex

echo "MacBook campaign completed. Review logs and environment metadata before release."
