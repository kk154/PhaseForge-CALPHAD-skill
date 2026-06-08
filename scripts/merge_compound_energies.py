#!/usr/bin/env python3
"""Merge compound energy CSV files, keeping the lowest energy per phase."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    best: dict[str, dict[str, str]] = {}
    fieldnames: list[str] | None = None
    for input_path in args.inputs:
        with Path(input_path).open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if fieldnames is None:
                fieldnames = list(reader.fieldnames or [])
            for row in reader:
                phase = row["phase"].upper()
                energy = float(row["energy_eV_atom"])
                if phase not in best or energy < float(best[phase]["energy_eV_atom"]):
                    row["phase"] = phase
                    row["source_table"] = input_path
                    best[phase] = row

    assert fieldnames is not None
    if "source_table" not in fieldnames:
        fieldnames.append("source_table")
    with Path(args.output).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(best.values(), key=lambda r: r["phase"]):
            writer.writerow(row)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

