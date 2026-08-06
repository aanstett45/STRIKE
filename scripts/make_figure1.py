#!/usr/bin/env python3
"""Generate the Figure 1 PGFPlots blocks from a full-SpGEMM CSV."""
import csv
import statistics
import sys
from pathlib import Path


def fail(message: str) -> "NoReturn":
    raise SystemExit(message)


if len(sys.argv) != 2:
    fail(f"Usage: {Path(sys.argv[0]).name} fullspgemm.csv")

path = Path(sys.argv[1])
if not path.is_file():
    fail(f"Input CSV not found: {path}")

with path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

required = {
    "kind",
    "param",
    "t_stateless_ms",
    "t_guarded_ms",
    "t_stateless_gallop_ms",
    "rewind_rate",
}
missing = required.difference(rows[0].keys() if rows else set())
if missing:
    fail(f"Missing CSV columns: {', '.join(sorted(missing))}")

betas: list[str] = []
for row in rows:
    if row["kind"] == "beta" and row["param"] not in betas:
        betas.append(row["param"])
if not betas:
    fail("No beta rows found")


def series(column: str) -> list[tuple[str, float]]:
    output = []
    for beta in betas:
        values = [
            float(row[column])
            for row in rows
            if row["kind"] == "beta" and row["param"] == beta
        ]
        if not values:
            fail(f"No values for beta={beta}, column={column}")
        output.append((beta, statistics.median(values)))
    return output


stateless = series("t_stateless_ms")
guarded = series("t_guarded_ms")
stateless_gallop = series("t_stateless_gallop_ms")
rewind_rate = series("rewind_rate")


def fmt(pairs: list[tuple[str, float]]) -> str:
    return " ".join(f"({beta},{value:.6g})" for beta, value in pairs)


print("\\addplot+[mark=o] coordinates")
print("  {" + fmt(rewind_rate) + "};")
print("\\addlegendentry{rewind rate}")
print("\\addplot+[mark=square*] coordinates")
print(
    "  {"
    + fmt(
        [
            (beta, s / g)
            for (beta, s), (_, g) in zip(stateless, guarded, strict=True)
        ]
    )
    + "};"
)
print("\\addlegendentry{stateless / persistent}")
print("\\addplot+[mark=triangle*] coordinates")
print(
    "  {"
    + fmt(
        [
            (beta, s / g)
            for (beta, s), (_, g) in zip(
                stateless, stateless_gallop, strict=True
            )
        ]
    )
    + "};"
)
print("\\addlegendentry{plain / galloping (stateless)}")
