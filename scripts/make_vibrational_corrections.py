#!/usr/bin/env python3
"""Convert vibrational free energies to formation free-energy corrections."""

from __future__ import annotations

import ast
import csv
from pathlib import Path


EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002


def counts_from_row(row: dict[str, str]) -> dict[str, float]:
    raw = ast.literal_eval(row["formula_counts"])
    return {str(k).upper().replace("0+", ""): float(v) for k, v in raw.items()}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fvib", required=True, help="vibrational_free_energies.csv")
    parser.add_argument("--references", required=True, help="CSV with columns element,phase")
    parser.add_argument("--output", required=True)
    parser.add_argument("--table-output", default="")
    args = parser.parse_args()

    with Path(args.fvib).open(newline="", encoding="utf-8") as f:
        fvib_rows = {row["phase"]: row for row in csv.DictReader(f)}

    reference_phase_by_element: dict[str, str] = {}
    with Path(args.references).open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            reference_phase_by_element[row["element"].strip().upper()] = row["phase"].strip()

    refs = {
        element: float(fvib_rows[phase]["fvib_eV_atom"])
        for element, phase in reference_phase_by_element.items()
    }

    rows = []
    for phase, row in fvib_rows.items():
        if phase in reference_phase_by_element.values():
            continue
        counts = counts_from_row(row)
        total = sum(counts.values())
        missing = [el for el in counts if el not in refs]
        if missing:
            raise ValueError(f"{phase} missing reference elements: {', '.join(sorted(missing))}")
        reference = sum(counts.get(el, 0.0) / total * refs[el] for el in counts)
        formation = float(row["fvib_eV_atom"]) - reference
        rows.append(
            {
                "phase": phase,
                "formula_counts": counts,
                "fvib_eV_atom": row["fvib_eV_atom"],
                "reference_fvib_eV_atom": reference,
                "fvib_form_eV_atom": formation,
                "fvib_form_J_mol_atom": formation * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                "imaginary_modes": row.get("imaginary_modes", ""),
            }
        )

    table_output = Path(args.table_output) if args.table_output else Path(args.output).with_name("vibrational_formation_energies.csv")
    with table_output.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys()) if rows else [
            "phase",
            "formula_counts",
            "fvib_eV_atom",
            "reference_fvib_eV_atom",
            "fvib_form_eV_atom",
            "fvib_form_J_mol_atom",
            "imaginary_modes",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with Path(args.output).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["phase", "correction_eV_atom", "correction_J_mol_atom"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "phase": row["phase"],
                    "correction_eV_atom": row["fvib_form_eV_atom"],
                    "correction_J_mol_atom": row["fvib_form_J_mol_atom"],
                }
            )
    print(f"Wrote {table_output}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
