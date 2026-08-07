#!/usr/bin/env python3
"""Combine independently generated per-sample diagnostic CSV files."""

import csv
from pathlib import Path
import sys


def combine(output: Path, inputs: list[Path]) -> None:
    rows = []
    fieldnames = None
    for path in inputs:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"Diagnostic columns differ in {path}")
            rows.extend(reader)

    if not fieldnames:
        raise ValueError("No diagnostic records were provided")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["sample"]))


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: finalize_diagnostics.py OUTPUT.csv SAMPLE.csv [SAMPLE.csv ...]"
        )
    combine(Path(sys.argv[1]), [Path(value) for value in sys.argv[2:]])


if __name__ == "__main__":
    main()
