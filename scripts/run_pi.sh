#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

for f in data/roadNet-CA.txt data/email-Enron.txt; do
  [[ -f "$f" ]] || { echo "Missing required dataset: $f" >&2; exit 1; }
done

export PI_GOVERNOR="${PI_GOVERNOR:-performance}"
: "${PI_POWER_NOTES:?Set PI_POWER_NOTES to the actual power-supply description}"
: "${PI_COOLING_NOTES:?Set PI_COOLING_NOTES to the actual cooling/enclosure description}"
export OMP_SCHEDULE=static
export OMP_PROC_BIND=close
export OMP_PLACES=cores

bash scripts/check_pi_governor.sh "$PI_GOVERNOR"
bash scripts/reset_platform_results.sh raspberry_pi
mkdir -p build results/raspberry_pi/threads

g++ -O3 -march=native -fopenmp -std=c++17 \
  src/finger_spgemm.cpp -o build/finger_spgemm

bash scripts/capture_pi_environment.sh results/raspberry_pi/environment.txt

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

./build/finger_spgemm triangles \
  results/raspberry_pi/triangles_real.csv \
  --runs 5 --threads 1 \
  data/roadNet-CA.txt data/email-Enron.txt \
  2>&1 | tee results/raspberry_pi/triangles_real.log

./build/finger_spgemm adaptive \
  results/raspberry_pi/adaptive.csv \
  --sample-pct 2 --threshold 0.05 --runs 5 \
  data/roadNet-CA.txt data/email-Enron.txt \
  gen:band:20000:16 gen:ba:20000:8 \
  2>&1 | tee results/raspberry_pi/adaptive.log

./build/finger_spgemm fullspgemm \
  results/raspberry_pi/fullspgemm_pi.csv 2048 3 3 --max-band 64 \
  2>&1 | tee results/raspberry_pi/fullspgemm_pi.log

bash scripts/check_pi_governor.sh "$PI_GOVERNOR"
bash scripts/capture_pi_environment.sh results/raspberry_pi/environment_after.txt
python3 scripts/validate_results.py --platform raspberry_pi --pi-governor "$PI_GOVERNOR"

echo "Raspberry Pi campaign completed and locally validated."
