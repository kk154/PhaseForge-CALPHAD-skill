#!/usr/bin/env python3
"""Relax CIF structures with ORB and save compound energies."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from pathlib import Path

from pymatgen.core import Composition, Structure
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation

from materialsframework.calculators import ORBCalculator


EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002


def formula_counts_from_name(path: Path) -> dict[str, float]:
    stem = path.stem.split("_")[-1]
    comp = Composition(stem)
    return {str(el).upper(): float(amount) for el, amount in comp.as_dict().items()}


def parse_components(text: str) -> list[str]:
    return [item.strip().upper() for item in text.split(",") if item.strip()]


def phase_name_from_formula(path: Path, components: list[str]) -> str:
    formula = path.stem.split("_")[-1]
    counts = {str(el).upper(): amount for el, amount in Composition(formula).as_dict().items()}
    parts = []
    order = components if components else sorted(counts)
    for el in order:
        amount = counts.get(el, 0)
        if amount:
            num = "" if abs(float(amount) - 1.0) < 1e-12 else f"{float(amount):g}".replace(".", "P")
            parts.append(f"{el}{num}")
    return "P_" + "".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif-dir", default="../..", help="Directory containing CIF files")
    parser.add_argument("--output-dir", default="compound_orb")
    parser.add_argument("--model", default="orb-v3-conservative-inf-omat")
    parser.add_argument("--device", default=os.environ.get("PHASEFORGE_ORB_DEVICE", "cuda"))
    parser.add_argument("--compile", action="store_true")
    parser.add_argument("--max-orderings", type=int, default=10)
    parser.add_argument("--components", default="", help="Optional phase-name element order, e.g. A,B,C")
    args = parser.parse_args()

    cif_dir = Path(args.cif_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    components = parse_components(args.components)

    calculator = ORBCalculator(
        model=args.model,
        fmax=0.001,
        verbose=True,
        device=args.device,
        compile=args.compile,
    )

    rows: list[dict[str, object]] = []

    for cif_path in sorted(cif_dir.glob("*.cif")):
        workdir = output_dir / cif_path.stem
        workdir.mkdir(exist_ok=True)
        result_json = workdir / "result.json"

        if result_json.exists():
            rows.append(json.loads(result_json.read_text(encoding="utf-8")))
            continue

        print(f"Relaxing {cif_path.name}")
        structure = Structure.from_file(str(cif_path))
        ordered_from_disordered = False
        candidate_structures = [structure]
        if not structure.is_ordered:
            ranked = OrderDisorderedStructureTransformation().apply_transformation(
                structure, return_ranked_list=args.max_orderings
            )
            ordered_from_disordered = True
            candidate_structures = [entry["structure"] for entry in ranked]

        best_result = None
        best_structure = None
        best_index = 0
        for idx, candidate in enumerate(candidate_structures):
            print(f"  candidate {idx + 1}/{len(candidate_structures)}")
            result = calculator.relax(candidate)
            energy_eV_cell = float(result["energy"])
            energy_eV_atom_candidate = energy_eV_cell / float(len(candidate))
            if best_result is None or energy_eV_atom_candidate < best_result["energy_eV_atom"]:
                best_result = {
                    "raw": result,
                    "energy_eV_cell": energy_eV_cell,
                    "energy_eV_atom": energy_eV_atom_candidate,
                }
                best_structure = result["final_structure"]
                best_index = idx

        assert best_result is not None and best_structure is not None
        structure = candidate_structures[best_index]
        final_structure = best_structure
        final_structure.to(workdir / "CONTCAR", fmt="poscar")
        final_structure.to(workdir / f"{cif_path.stem}_relaxed.cif", fmt="cif")

        formula_counts = formula_counts_from_name(cif_path)
        formula_atoms = sum(formula_counts.values())
        structure_atoms = float(len(structure))
        energy_eV_cell = best_result["energy_eV_cell"]
        energy_eV_atom = best_result["energy_eV_atom"]

        row = {
            "source_cif": str(cif_path),
            "phase": phase_name_from_formula(cif_path, components),
            "display_formula": cif_path.stem,
            "parsed_formula": structure.composition.reduced_formula,
            "ordered_from_disordered": ordered_from_disordered,
            "selected_ordering_index": best_index,
            "num_orderings_tested": len(candidate_structures),
            "formula_counts": formula_counts,
            "formula_atoms": formula_atoms,
            "structure_atoms": structure_atoms,
            "energy_eV_cell": energy_eV_cell,
            "energy_eV_atom": energy_eV_atom,
            "energy_J_mol_atom": energy_eV_atom * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
            "model": args.model,
        }
        result_json.write_text(json.dumps(row, indent=2), encoding="utf-8")
        rows.append(row)

    csv_path = output_dir / "compound_energies.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "source_cif",
                "phase",
                "display_formula",
                "parsed_formula",
                "ordered_from_disordered",
                "selected_ordering_index",
                "num_orderings_tested",
                "formula_counts",
                "formula_atoms",
                "structure_atoms",
                "energy_eV_cell",
                "energy_eV_atom",
                "energy_J_mol_atom",
                "model",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
