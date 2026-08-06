[![DOI](https://zenodo.org/badge/1113783061.svg)](https://doi.org/10.5281/zenodo.21819709)

# Cross-row Pointer Persistence in Merge-based Sparse Matrix Products

Reproducibility package for the article **“Cross-row pointer persistence in merge-based sparse matrix products.”**

This tree is the corrected **v1.0.1 rerun package**. The benchmark outputs from the earlier archive were removed because they predated the portability fixes and the Raspberry Pi environment was inconsistent. Run both campaigns again before creating the final v1.0.1 tag and Zenodo release.

## Important v1.0.1 corrections

- Synthetic inputs no longer depend on `std::shuffle` or `std::uniform_int_distribution`; Mac/libc++ and Pi/libstdc++ now generate identical inputs.
- RCM ordering breaks degree ties by vertex ID and is reproducible across standard libraries.
- Every generated graph and full-SpGEMM input has a stable 64-bit structure fingerprint recorded in CSV files and logs.
- The duplicated Full-SpGEMM CSV header fields were removed; each data row now has exactly the same number of columns as the header.
- Invalid run counts, thread counts, sampling percentages, thresholds, matrix sizes, and band widths are rejected.
- Timings use `std::chrono::steady_clock`; even-sized timing series use the mathematical median.
- Adaptive sampling is deterministically spread over the full row range instead of using only the first rows.
- Mac environment capture removes serial/UUID/user/hostname fields. Pi capture records the actual model, memory, bitness, governors, temperature, and throttling state.
- The Pi campaign refuses to run unless the requested governor is active.
- Complete scripts now generate `triangles_real.csv`, `adaptive.csv`, and `fullspgemm_pi.csv`.

See [`CHANGELOG.md`](CHANGELOG.md) for the complete change list.

## Hardware for the new measurements

### MacBook

- MacBook Pro with Apple M4 Pro
- 24 GB memory
- macOS environment captured at run time
- OpenMP schedule: `static`, except the explicit scheduling comparison

### Raspberry Pi

The archived hardware evidence identifies the machine as:

- Raspberry Pi 3 Model B Plus Rev 1.3
- Broadcom BCM2837B0 / Cortex-A53, four in-order cores
- maximum CPU frequency 1.4 GHz
- 1 GB nominal memory, approximately 921 MiB visible to the operating system
- `armv7l`, 32-bit userland

The corrected campaign defaults to the `performance` governor and verifies it before and after the run. To intentionally use another governor, set `PI_GOVERNOR` consistently and report that actual setting in the article.

## Repository layout

```text
src/finger_spgemm.cpp          benchmark and verification program
scripts/run_macbook.sh         complete Mac campaign
scripts/run_pi.sh              complete Raspberry Pi campaign
scripts/set_pi_governor.sh     set and verify the Pi governor
scripts/validate_results.py    per-platform and cross-platform validation
scripts/prepare_release.sh     final checks and manifest regeneration
scripts/package_release.sh     validated v1.0.1 ZIP creation
data/                          SNAP input datasets and provenance
results/                       regenerated raw CSV, logs, and environments
generated/                     regenerated article fragments
CHANGELOG.md                   detailed corrections
RERUN_GUIDE.md                 exact rerun workflow
```

## Prerequisites

### macOS

```bash
xcode-select --install
brew install libomp python
```

### Raspberry Pi OS / Linux

```bash
sudo apt update
sudo apt install -y build-essential python3 curl gzip
```

## Datasets

The two SNAP files are already included in this package. Verify them from the repository root:

```bash
sha256sum -c data/checksums.sha256      # Linux
shasum -a 256 -c data/checksums.sha256  # macOS
```

They can also be downloaded again with:

```bash
bash scripts/download_snap_data.sh
```

## Recommended rerun workflow

Use the **same repository tree** for both machines so the Mac results travel with the source to the Pi and return with the Pi results.

### 1. Run the Mac campaign

```bash
bash scripts/run_macbook.sh
```

It produces:

- verification log;
- thread sweep for 1, 2, 4, 8, and 12 threads;
- sequential real-graph campaign;
- OpenMP scheduling comparison;
- adaptive campaign;
- main and extended Full-SpGEMM campaigns;
- regenerated Figure 1 fragment.

### 2. Move the same tree to the Raspberry Pi

Set the governor first:

```bash
sudo bash scripts/set_pi_governor.sh performance
```

Record concrete power and cooling conditions when starting the run:

```bash
PI_POWER_NOTES="official 5.1 V power supply" \
PI_COOLING_NOTES="passive heatsink, open air" \
bash scripts/run_pi.sh
```

Replace those descriptions with the actual setup. The Pi campaign produces:

- environment captures before and after the measurements;
- verification log;
- thread sweep for 1, 2, and 4 threads;
- sequential real-graph campaign;
- adaptive campaign;
- `fullspgemm_pi.csv` with the beta sweep and band widths through 64.

### 3. Validate the combined result set

After bringing the completed repository back to the main machine:

```bash
python3 scripts/validate_results.py \
  --platform all \
  --compare-platforms \
  --pi-governor performance
```

This checks required files, CSV shape, correctness flags, verification logs, privacy fields, Pi identity/governor, and identical Mac/Pi input fingerprints.

### 4. Prepare the final release

```bash
bash scripts/prepare_release.sh
```

Then review the article against the regenerated raw values. When the article, repository, and metadata agree:

```bash
bash scripts/package_release.sh
```

The package script writes `STRIKE-v1.0.1.zip` one directory above the repository.

## Compilation commands

### MacBook

```bash
mkdir -p build
clang++ -O3 -march=native \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -std=c++17 src/finger_spgemm.cpp \
  -o build/finger_spgemm
```

### Raspberry Pi

```bash
mkdir -p build
g++ -O3 -march=native -fopenmp -std=c++17 \
  src/finger_spgemm.cpp \
  -o build/finger_spgemm
```

## Verification values for the corrected portable generators

The exact non-guarded results below are expected on both standard libraries after compiling this source:

| Configuration | Difference count | Survival rate | Input fingerprints |
|---|---:|---:|---|
| Random `n=512`, `d=0.02`, `beta=0.6` | 45,391 | 0.759233% | `6355167ee4b9ce8f` / `422484e9fc6f82ad` |
| Band `n=1024`, `w=16` | 30,255 | 51.2772% | `b3f888ce0d70b23d` / same |
| Random `n=2048`, `d=0.01`, `beta=0.6` | 736,847 | 0.104397% | `43f00c81ebdaf5d0` / `8f2e7e02ca145649` |

The stateless, guarded, and guarded-plus-galloping kernels must match the Gustavson reference. Every triangle cross-check must end with `[OK]`.

Because the corrected `n=2048` value is slightly above 0.1%, the article should retain approximate wording such as “roughly one reference nonzero in a thousand,” not a strict “under 0.1%” claim.

## Full-SpGEMM options

```text
finger_spgemm fullspgemm OUT.csv [n] [trials] [runs] [--max-band W]
```

The Mac scripts use `--max-band 512`. The Pi script uses `--max-band 64`, so the run finishes cleanly without manual interruption or a partial CSV.

## Raw-output policy

Keep CSV and `.log` files unmodified. Logs preserve every `[spread]` line. Do not replace raw outputs with rounded tables or spreadsheet exports. If a campaign is rerun, the previous platform directory is moved to `.benchmark_backups/`, which must be removed or stored elsewhere before packaging the release.

## License

The code and repository documentation are released under the MIT License. Dataset files remain subject to the terms and attribution requirements of their original providers.
