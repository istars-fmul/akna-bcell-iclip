#!/usr/bin/env python3
"""Convert a Flexbar length distribution into iCLIPro group arguments."""

from pathlib import Path
import sys


def parse_length_distribution(path: Path) -> tuple[str, str]:
    lengths = set()
    with path.open(encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            fields = line.split()
            if fields:
                lengths.add(int(fields[0]))

    if len(lengths) < 2:
        raise ValueError(f"Expected at least two read lengths in {path}")

    ordered = sorted(lengths)
    groups = ",".join(
        [f"L{length}:{length}" for length in ordered[:-1]]
        + [f"R:{ordered[-1]}"]
    )
    comparisons = ",".join(f"L{length}-R" for length in ordered[:-1])
    return groups, comparisons


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: parse_length_distribution.py LENGTHDIST")
    print(*parse_length_distribution(Path(sys.argv[1])))


if __name__ == "__main__":
    main()
