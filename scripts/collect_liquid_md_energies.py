#!/usr/bin/env python3
"""Collect LIQUID LAMMPS/MLIP MD summaries into one CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def composition_text(composition: dict[str, float]) -> str:
    return ";".join(f"{el}:{composition[el]:.12g}" for el in sorted(composition))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="LIQUID")
    parser.add_argument("--output", default="liquid_md_energies.csv")
    args = parser.parse_args()

    root = Path(args.root)
    summaries = sorted(root.rglob("liquid_md_summary.json"))
    if not summaries:
        raise FileNotFoundError(f"No liquid_md_summary.json files found under {root}")

    rows: list[dict[str, object]] = []
    for path in summaries:
        data = json.loads(path.read_text(encoding="utf-8"))
        composition = {str(k).upper(): float(v) for k, v in data["composition"].items()}
        rows.append(
            {
                "phase": data.get("phase", "LIQUID"),
                "composition": composition_text(composition),
                "composition_json": json.dumps(composition, sort_keys=True),
                "formula_counts": json.dumps(data.get("formula_counts", {}), sort_keys=True),
                "temperature_K": data["temperature_K"],
                "pe_eV_atom": data["pe_eV_atom"],
                "pv_eV_atom": data["pv_eV_atom"],
                "enthalpy_eV_atom": data["enthalpy_eV_atom"],
                "press_bar": data.get("press_bar", ""),
                "volume_A3": data.get("volume_A3", ""),
                "source_dir": str(path.parent),
            }
        )

    fieldnames = [
        "phase",
        "composition",
        "composition_json",
        "formula_counts",
        "temperature_K",
        "pe_eV_atom",
        "pv_eV_atom",
        "enthalpy_eV_atom",
        "press_bar",
        "volume_A3",
        "source_dir",
    ]
    with Path(args.output).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output} with {len(rows)} rows")


if __name__ == "__main__":
    main()
