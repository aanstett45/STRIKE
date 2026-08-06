# Changelog

## v1.0.1 — corrected rerun package

### Release and archive

- Updated package metadata from v1.0.0 to v1.0.1.
- Removed stale benchmark outputs and the stale generated Figure 1 fragment. They must be regenerated from the corrected source.
- Added an exact rerun guide, result validator, release preparation script, packaging script, and reproducible manifest generator.
- Added automatic backup of a platform’s previous outputs to `.benchmark_backups/` before a new complete campaign.

### Raspberry Pi identity and measurement state

- Corrected all repository descriptions from Raspberry Pi 4 / 4 GB to Raspberry Pi 3 Model B Plus Rev 1.3 / 1 GB nominal / approximately 921 MiB visible / armv7l 32-bit / 1.4 GHz.
- Removed hard-coded false hardware values from the Pi environment capture; model and memory are now read from the machine.
- Added `set_pi_governor.sh` and `check_pi_governor.sh`.
- The complete Pi campaign now refuses to start or finish with a governor different from `PI_GOVERNOR` (`performance` by default).
- Added environment captures before and after the Pi campaign to record governor, temperature, and throttling state.
- Replaced release placeholders with optional environment variables for actual power and cooling notes.

### Privacy

- Mac environment capture now removes serial number, hardware UUID, provisioning UDID, computer name, user name, and Activation Lock fields.
- Replaced `uname -a` with `uname -smr` so hostnames are not recorded.
- Removed old environment captures containing personal machine and account identifiers.

### Portable and deterministic inputs

- Replaced `std::uniform_int_distribution` with a portable rejection-sampling helper based only on specified `mt19937_64` outputs.
- Replaced `std::shuffle` with a portable Fisher–Yates implementation.
- Added deterministic degree/vertex tie-breaking to both RCM sorts.
- Added stable structure fingerprints to verification output and benchmark CSV files.
- Added cross-platform validation that rejects Mac/Pi campaigns using different generated or loaded graph structures.

### Benchmark correctness and robustness

- Fixed the Full-SpGEMM CSV header, which previously duplicated two timing names and had more columns than each data row.
- Added output-file creation and write-error checks.
- Added argument validation for sizes, run counts, trial counts, thread counts, sampling percentages, thresholds, and band widths.
- Added explicit errors for unknown or incomplete command-line options.
- Changed benchmark timing from `high_resolution_clock` to monotonic `steady_clock`.
- Fixed the median for even-sized timing series.
- Changed adaptive sampling from the first rows to deterministic rows distributed across the entire graph.
- Added `--max-band` so the Pi campaign can stop cleanly at band width 64.
- Removed compiler warnings caused by unused formal parameters and the line-continuation character in a comment.

### Campaign coverage

The Mac complete campaign now generates:

- `results/macbook/triangles_real.csv` and `.log`;
- `results/macbook/adaptive.csv` and `.log`;
- the existing thread, scheduling, main Full-SpGEMM, and extended-band outputs.

The Raspberry Pi complete campaign now generates:

- `results/raspberry_pi/triangles_real.csv` and `.log`;
- `results/raspberry_pi/adaptive.csv` and `.log`;
- `results/raspberry_pi/fullspgemm_pi.csv` and `.log`;
- environment captures before and after the campaign.

### Figure generation

- Replaced the custom incorrect even-sized median with Python’s `statistics.median`.
- Added argument, file, header, and data validation to `make_figure1.py`.
- Increased output precision so generated coordinates are not unnecessarily rounded to three decimals.

### New portable verification baseline

The corrected source produces the following exact non-guarded Full-SpGEMM values on the tested Linux/libstdc++ build and, by construction, identical input fingerprints across conforming implementations:

- random 512: 45,391 differences, 0.759233% survival;
- band 1024: 30,255 differences, 51.2772% survival;
- random 2048: 736,847 differences, 0.104397% survival.

Both machines must be rerun. Values from the pre-v1.0.1 CSV and logs must not be mixed with the corrected campaign.
