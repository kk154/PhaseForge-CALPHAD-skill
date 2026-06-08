#!/usr/bin/env python3
"""Prepare phonon target and reference CSV files from relaxed structures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def truthy(text: str) -> bool:
    return text.strip().lower() in {"1", "true", "yes", "y"}


def read_reference_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    required = {"element", "phase", "path"}
    missing = required.difference(rows[0].keys() if rows else [])
    if missing:
        raise ValueError(f"Reference CSV missing columns: {', '.join(sorted(missing))}")
    return rows


def compound_structure_path(row: dict[str, str], compound_dir: Path) -> Path:
    if row.get("relaxed_structure"):
        return Path(row["relaxed_structure"])
    display = row.get("display_formula", "").strip()
    if not display:
        source = Path(row["source_cif"])
        display = source.stem
    return compound_dir / display / "CONTCAR"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compound-energies", default="compound_orb/compound_energies.csv")
    parser.add_argument("--compound-dir", default="compound_orb")
    parser.add_argument(
        "--reference-structures",
        required=True,
        help="CSV with columns element,phase,path for elemental reference structures",
    )
    parser.add_argument("--targets-output", default="phonon_targets.csv")
    parser.add_argument("--references-output", default="phonon_references.csv")
    parser.add_argument(
        "--phases",
        nargs="*",
        help="Optional compound phase subset. Element references are always included.",
    )
    parser.add_argument(
        "--include-reference-targets",
        default="true",
        help="Set false only when reference phases already exist in the fvib CSV.",
    )
    args = parser.parse_args()

    requested = {phase.upper() for phase in args.phases} if args.phases else None
    compound_dir = Path(args.compound_dir)
    target_rows: list[dict[str, str]] = []
    reference_rows = read_reference_rows(Path(args.reference_structures))

    if truthy(args.include_reference_targets):
        for row in reference_rows:
            path = Path(row["path"])
            if not path.exists():
                raise FileNotFoundError(f"Reference {row['phase']}: {path}")
            target_rows.append({"phase": row["phase"].strip(), "path": str(path), "is_reference": "true"})

    with Path(args.compound_energies).open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            phase = row["phase"].strip()
            if requested and phase.upper() not in requested:
                continue
            path = compound_structure_path(row, compound_dir)
            if not path.exists():
                raise FileNotFoundError(f"{phase}: relaxed structure not found at {path}")
            target_rows.append({"phase": phase, "path": str(path), "is_reference": "false"})

    with Path(args.targets_output).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["phase", "path", "is_reference"])
        writer.writeheader()
        writer.writerows(target_rows)

    with Path(args.references_output).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["element", "phase"])
        writer.writeheader()
        for row in reference_rows:
            writer.writerow({"element": row["element"].strip().upper(), "phase": row["phase"].strip()})

    print(f"Wrote {args.targets_output}")
    print(f"Wrote {args.references_output}")


if __name__ == "__main__":
    main()
