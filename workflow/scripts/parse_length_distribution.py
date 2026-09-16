#!/usr/bin/env python3
"""Convert a Flexbar length distribution into iCLIPro group arguments.

Input: a whitespace-delimited .lengthdist file with a header and read length in
the first column. Counts are not used to define groups; each listed length is
included once. At least two lengths are needed for a reference comparison.

Output on stdout: two space-separated strings for iCLIPro's -g and -p options.
For lengths 15, 16 and 17, these are L15:15,L16:16,R:17 and L15-R,L16-R.
The workflow captures this line in a temporary file for the iCLIPro rule.
"""

from pathlib import Path
import sys


def parse_length_distribution(path: Path) -> tuple[str, str]:
    # Collect unique read lengths from the Flexbar table, skipping its header.
    lengths = set()
    with path.open(encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            fields = line.split()
            if fields:
                lengths.add(int(fields[0]))

    if len(lengths) < 2:
        raise ValueError(f"Expected at least two read lengths in {path}")

    # Use the longest listed length as reference R; label the others L<length>.
    ordered = sorted(lengths)
    groups = ",".join(
        [f"L{length}:{length}" for length in ordered[:-1]]
        + [f"R:{ordered[-1]}"]
    )
    # Compare every non-reference length group with the common reference.
    comparisons = ",".join(f"L{length}-R" for length in ordered[:-1])
    return groups, comparisons


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: parse_length_distribution.py LENGTHDIST")
    print(*parse_length_distribution(Path(sys.argv[1])))


if __name__ == "__main__":
    main()
