# Raw benchmark results

This directory is intentionally empty in the corrected rerun package, apart from `.gitkeep` files. Do not restore pre-v1.0.1 measurements: the synthetic generators, RCM ordering, CSV schema, and Pi measurement controls changed.

The complete scripts generate the following files.

## MacBook

- `environment.txt`
- `verify.log`
- `triangles_real.csv` and `.log`
- `adaptive.csv` and `.log`
- `fullspgemm.csv` and `.log`
- `full_wide.csv` and `.log`
- `threads/tri_T{1,2,4,8,12}.csv` and `.log`
- `scheduling/{static,dynamic64,guided}.csv` and `.log`

## Raspberry Pi 3 B+

- `environment.txt`
- `environment_after.txt`
- `verify.log`
- `triangles_real.csv` and `.log`
- `adaptive.csv` and `.log`
- `fullspgemm_pi.csv` and `.log`
- `threads/tri_T{1,2,4}.csv` and `.log`

Keep raw CSV and logs unmodified. Run `python3 scripts/validate_results.py --platform all --compare-platforms` after both campaigns.
