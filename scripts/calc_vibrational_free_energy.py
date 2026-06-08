#!/usr/bin/env python3
"""Calculate ORB/phonopy vibrational free energies for a target list."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np
from ase import Atoms
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure

from materialsframework.calculators import ORBCalculator


EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002


def phonopy_to_ase(atoms: PhonopyAtoms) -> Atoms:
    return Atoms(
        symbols=atoms.symbols,
        cell=atoms.cell,
        scaled_positions=atoms.scaled_positions,
        pbc=True,
    )


def formula_counts(path: Path) -> dict[str, float]:
    structure = Structure.from_file(str(path))
    return {str(el).upper().replace("0+", ""): float(amount) for el, amount in structure.composition.as_dict().items()}


def reduced_counts(counts: dict[str, float]) -> dict[str, float]:
    values = [int(round(v)) for v in counts.values()]
    gcd = values[0]
    for value in values[1:]:
        gcd = math.gcd(gcd, value)
    if gcd <= 0:
        return counts
    return {el: value / gcd for el, value in counts.items()}


def thermal_free_energy_at_temperature(thermal_properties: dict, temperature: float) -> float:
    temperatures = np.asarray(thermal_properties["temperatures"], dtype=float)
    free_energy = np.asarray(thermal_properties["free_energy"], dtype=float)
    return float(np.interp(temperature, temperatures, free_energy))


def parse_triplet(text: str) -> list[int]:
    values = [int(x) for x in text.split(",")]
    if len(values) != 3:
        raise ValueError("Expected three comma-separated integers")
    return values


def calculate_fvib(
    structure_path: Path,
    output_dir: Path,
    calculator: ORBCalculator,
    supercell_matrix: list[list[int]],
    displacement_distance: float,
    mesh: list[int],
    temperature: float,
    tmax: float,
    tstep: float,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    result_json = output_dir / "phonon_result.json"
    if result_json.exists():
        return json.loads(result_json.read_text(encoding="utf-8"))

    unitcell, _ = read_crystal_structure(str(structure_path), interface_mode="vasp")
    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, primitive_matrix="P")
    phonon.generate_displacements(distance=displacement_distance)
    supercells = phonon.supercells_with_displacements

    forces = []
    for idx, supercell in enumerate(supercells, start=1):
        displaced_dir = output_dir / f"disp-{idx:03d}"
        displaced_dir.mkdir(exist_ok=True)
        force_path = displaced_dir / "forces.npy"
        if force_path.exists():
            force = np.load(force_path)
        else:
            atoms = phonopy_to_ase(supercell)
            result = calculator.calculate(atoms)
            force = np.asarray(result["forces"], dtype=float)
            np.save(force_path, force)
        forces.append(force)

    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.run_mesh(mesh, is_mesh_symmetry=False)
    mesh_dict = phonon.get_mesh_dict()
    frequencies = np.asarray(mesh_dict["frequencies"], dtype=float)
    min_frequency_thz = float(np.nanmin(frequencies))
    imaginary_modes = int(np.count_nonzero(frequencies < -1e-5))
    total_modes = int(frequencies.size)

    phonon.run_thermal_properties(t_step=tstep, t_max=tmax, t_min=0)
    thermal = phonon.get_thermal_properties_dict()
    phonon.save(output_dir / "phonopy_params.yaml")

    counts = formula_counts(structure_path)
    natoms_unitcell = float(sum(counts.values()))
    fvib_kj_mol_cell = thermal_free_energy_at_temperature(thermal, temperature)
    fvib_eV_cell = fvib_kj_mol_cell / 96.48533212331002
    fvib_eV_atom = fvib_eV_cell / natoms_unitcell

    result = {
        "structure": str(structure_path),
        "formula_counts": counts,
        "reduced_counts": reduced_counts(counts),
        "natoms_unitcell": natoms_unitcell,
        "num_displacements": len(supercells),
        "supercell_matrix": supercell_matrix,
        "mesh": mesh,
        "min_frequency_THz": min_frequency_thz,
        "imaginary_modes": imaginary_modes,
        "total_mesh_modes": total_modes,
        "temperature_K": temperature,
        "fvib_eV_atom": fvib_eV_atom,
        "fvib_J_mol_atom": fvib_eV_atom * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
    }
    result_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def load_targets(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    required = {"phase", "path"}
    missing = required.difference(rows[0].keys() if rows else [])
    if missing:
        raise ValueError(f"Targets CSV missing columns: {', '.join(sorted(missing))}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True, help="CSV with columns phase,path and optional is_reference")
    parser.add_argument("--output-dir", default="phonon_free_energy")
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--tmax", type=float, default=1400)
    parser.add_argument("--tstep", type=float, default=10)
    parser.add_argument("--displacement", type=float, default=0.01)
    parser.add_argument(
        "--supercell",
        default="2,2,2",
        help="Diagonal supercell. Default 2,2,2 is the production setting; use 1,1,1 only for fast screening.",
    )
    parser.add_argument("--mesh", default="12,12,12", help="Phonon q mesh")
    parser.add_argument("--phases", nargs="*", help="Optional subset of phase labels to calculate")
    parser.add_argument("--model", default="orb-v3-conservative-inf-omat")
    parser.add_argument("--device", default=os.environ.get("PHASEFORGE_ORB_DEVICE", "cuda"))
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    supercell_diag = parse_triplet(args.supercell)
    mesh = parse_triplet(args.mesh)
    supercell_matrix = [[supercell_diag[0], 0, 0], [0, supercell_diag[1], 0], [0, 0, supercell_diag[2]]]

    targets = load_targets(Path(args.targets))
    if args.phases:
        requested = {phase.upper() for phase in args.phases}
        targets = [row for row in targets if row["phase"].upper() in requested]
        missing = requested.difference({row["phase"].upper() for row in targets})
        if missing:
            raise ValueError(f"Unknown phases requested: {', '.join(sorted(missing))}")

    calculator = ORBCalculator(
        model=args.model,
        fmax=0.001,
        verbose=True,
        device=args.device,
        compile=False,
    )

    rows: list[dict[str, object]] = []
    for target in targets:
        phase = target["phase"]
        path = Path(target["path"])
        if not path.exists():
            raise FileNotFoundError(f"{phase}: {path}")
        print(f"Calculating phonons for {phase}: {path}")
        row = calculate_fvib(
            path,
            output_dir / phase,
            calculator,
            supercell_matrix,
            args.displacement,
            mesh,
            args.temperature,
            args.tmax,
            args.tstep,
        )
        row["phase"] = phase
        row["is_reference"] = target.get("is_reference", "").strip().lower() in {"1", "true", "yes", "y"}
        rows.append(row)

    fieldnames = [
        "phase",
        "is_reference",
        "structure",
        "formula_counts",
        "reduced_counts",
        "natoms_unitcell",
        "num_displacements",
        "supercell_matrix",
        "mesh",
        "min_frequency_THz",
        "imaginary_modes",
        "total_mesh_modes",
        "temperature_K",
        "fvib_eV_atom",
        "fvib_J_mol_atom",
    ]
    csv_path = output_dir / "vibrational_free_energies.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fieldnames} for row in rows)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
