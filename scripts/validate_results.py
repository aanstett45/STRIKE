#!/usr/bin/env python3
"""Validate benchmark outputs and cross-platform input identity."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MAC_EXPECTED = [
    "environment.txt",
    "verify.log",
    "triangles_real.csv",
    "triangles_real.log",
    "adaptive.csv",
    "adaptive.log",
    "fullspgemm.csv",
    "fullspgemm.log",
    "full_wide.csv",
    "full_wide.log",
    *[f"threads/tri_T{k}.{ext}" for k in (1, 2, 4, 8, 12) for ext in ("csv", "log")],
    *[f"scheduling/{name}.{ext}" for name in ("static", "dynamic64", "guided") for ext in ("csv", "log")],
]

PI_EXPECTED = [
    "environment.txt",
    "environment_after.txt",
    "verify.log",
    "triangles_real.csv",
    "triangles_real.log",
    "adaptive.csv",
    "adaptive.log",
    "fullspgemm_pi.csv",
    "fullspgemm_pi.log",
    *[f"threads/tri_T{k}.{ext}" for k in (1, 2, 4) for ext in ("csv", "log")],
]

SENSITIVE_PATTERNS = [
    re.compile(r"Serial Number:", re.I),
    re.compile(r"Hardware UUID:", re.I),
    re.compile(r"Provisioning UDID:", re.I),
    re.compile(r"Computer Name:", re.I),
    re.compile(r"User Name:", re.I),
]


class ValidationError(RuntimeError):
    pass


def require_files(platform: str) -> None:
    base = ROOT / "results" / platform
    expected = MAC_EXPECTED if platform == "macbook" else PI_EXPECTED
    missing = [str(base / rel) for rel in expected if not (base / rel).is_file()]
    empty = [str(base / rel) for rel in expected if (base / rel).is_file() and (base / rel).stat().st_size == 0]
    if missing:
        raise ValidationError("Missing files:\n  " + "\n  ".join(missing))
    if empty:
        raise ValidationError("Empty files:\n  " + "\n  ".join(empty))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValidationError(f"CSV has no header: {path}")
        rows = list(reader)
    if not rows:
        raise ValidationError(f"CSV has no data rows: {path}")
    for line_no, row in enumerate(rows, start=2):
        if None in row:
            raise ValidationError(f"Too many fields in {path}:{line_no}")
        if any(value is None for value in row.values()):
            raise ValidationError(f"Too few fields in {path}:{line_no}")
    return rows


def validate_csvs(platform: str) -> None:
    base = ROOT / "results" / platform
    for path in sorted(base.rglob("*.csv")):
        rows = read_csv(path)
        fields = set(rows[0])
        if "counts_all_equal" in fields:
            bad = [i for i, row in enumerate(rows, start=2) if row["counts_all_equal"] not in {"1", "true", "True"}]
            if bad:
                raise ValidationError(f"Triangle correctness failure in {path}, rows {bad}")
        if "correct_guarded_gallop" in fields:
            bad = [i for i, row in enumerate(rows, start=2) if row["correct_guarded_gallop"] not in {"1", "true", "True"}]
            if bad:
                raise ValidationError(f"SpGEMM correctness failure in {path}, rows {bad}")
        for fp_name in ("graph_fingerprint", "fingerprint_A", "fingerprint_B"):
            if fp_name in fields:
                for line_no, row in enumerate(rows, start=2):
                    if not re.fullmatch(r"[0-9a-f]{16}", row[fp_name]):
                        raise ValidationError(f"Invalid {fp_name} in {path}:{line_no}")


def validate_logs(platform: str) -> None:
    base = ROOT / "results" / platform
    verify = (base / "verify.log").read_text(encoding="utf-8", errors="replace")
    if "[ECHEC]" in verify or "garde+gallop=0" in verify:
        raise ValidationError(f"Verification failed in {base / 'verify.log'}")
    if verify.count("[OK]") < 3:
        raise ValidationError(f"Expected at least three triangle [OK] records in {base / 'verify.log'}")
    for path in base.rglob("*.log"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "ERREUR:" in text or "counts_ok=0" in text or "correct=0" in text:
            raise ValidationError(f"Failure marker found in {path}")


def validate_environment(platform: str, pi_governor: str) -> None:
    base = ROOT / "results" / platform
    for path in base.glob("environment*.txt"):
        text = path.read_text(encoding="utf-8", errors="replace")
        for pattern in SENSITIVE_PATTERNS:
            if pattern.search(text):
                raise ValidationError(f"Sensitive field found in {path}: {pattern.pattern}")
        if "EDIT THIS LINE" in text:
            raise ValidationError(f"Release placeholder found in {path}")
    if platform == "raspberry_pi":
        for path in (base / "environment.txt", base / "environment_after.txt"):
            text = path.read_text(encoding="utf-8", errors="replace")
            if "Raspberry Pi 3 Model B Plus" not in text:
                raise ValidationError(f"Unexpected Raspberry Pi model in {path}")
            governors = re.findall(r"scaling_governor:\s*(\S+)", text)
            if not governors or any(value != pi_governor for value in governors):
                raise ValidationError(f"Governor mismatch in {path}: {governors}")
            if "Power supply: not provided" in text or "Cooling: not provided" in text:
                raise ValidationError(f"Missing concrete Pi power/cooling notes in {path}")


def triangle_identity(platform: str) -> dict[tuple[str, str], tuple[str, str, str, str]]:
    base = ROOT / "results" / platform
    mapping: dict[tuple[str, str], tuple[str, str, str, str]] = {}
    for path in base.rglob("*.csv"):
        rows = read_csv(path)
        if "graph_fingerprint" not in rows[0]:
            continue
        for row in rows:
            key = (row["graph"], row.get("ordering", "adaptive"))
            value = (
                row["graph_fingerprint"],
                row.get("n", ""),
                row.get("nnz", ""),
                row.get("triangles", ""),
            )
            old = mapping.get(key)
            if old is not None and old != value:
                raise ValidationError(f"Inconsistent graph identity within {platform}: {key}: {old} vs {value}")
            mapping[key] = value
    return mapping


def full_identity(platform: str) -> dict[tuple[str, str, str], tuple[str, str]]:
    name = "fullspgemm.csv" if platform == "macbook" else "fullspgemm_pi.csv"
    rows = read_csv(ROOT / "results" / platform / name)
    mapping: dict[tuple[str, str, str], tuple[str, str]] = {}
    for row in rows:
        key = (row["kind"], row["param"], row["n"])
        value = (row["fingerprint_A"], row["fingerprint_B"])
        old = mapping.get(key)
        if old is not None and old != value:
            # beta trials intentionally have different seeds; include each pair as a multiset later.
            continue
        mapping[key] = value
    return mapping


def full_identity_multiset(platform: str) -> dict[tuple[str, str, str], list[tuple[str, str]]]:
    name = "fullspgemm.csv" if platform == "macbook" else "fullspgemm_pi.csv"
    rows = read_csv(ROOT / "results" / platform / name)
    mapping: dict[tuple[str, str, str], list[tuple[str, str]]] = {}
    for row in rows:
        key = (row["kind"], row["param"], row["n"])
        mapping.setdefault(key, []).append((row["fingerprint_A"], row["fingerprint_B"]))
    return mapping


def compare_platforms() -> None:
    mac_tri = triangle_identity("macbook")
    pi_tri = triangle_identity("raspberry_pi")
    common = sorted(set(mac_tri).intersection(pi_tri))
    if not common:
        raise ValidationError("No common triangle inputs found between platforms")
    for key in common:
        if mac_tri[key] != pi_tri[key]:
            raise ValidationError(f"Mac/Pi triangle input mismatch for {key}: {mac_tri[key]} vs {pi_tri[key]}")

    mac_full = full_identity_multiset("macbook")
    pi_full = full_identity_multiset("raspberry_pi")
    common_full = sorted(set(mac_full).intersection(pi_full))
    if not common_full:
        raise ValidationError("No common full-SpGEMM inputs found between platforms")
    for key in common_full:
        if mac_full[key] != pi_full[key]:
            raise ValidationError(f"Mac/Pi full-SpGEMM input mismatch for {key}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--platform", choices=("macbook", "raspberry_pi", "all"), default="all")
    parser.add_argument("--compare-platforms", action="store_true")
    parser.add_argument("--pi-governor", default="performance")
    args = parser.parse_args()

    platforms = ("macbook", "raspberry_pi") if args.platform == "all" else (args.platform,)
    try:
        for platform in platforms:
            require_files(platform)
            validate_csvs(platform)
            validate_logs(platform)
            validate_environment(platform, args.pi_governor)
            print(f"OK: {platform} result set")
        if args.compare_platforms:
            compare_platforms()
            print("OK: Mac/Pi input fingerprints match")
    except ValidationError as exc:
        print(f"VALIDATION ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
