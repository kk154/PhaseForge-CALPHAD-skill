#!/usr/bin/env python3
"""Append stoichiometric compound phases to a PhaseForge TDB."""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
from pathlib import Path

from pymatgen.core import Structure


R_J_MOL_K = 8.31446261815324
EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002


def parse_components(text: str) -> list[str]:
    components = [item.strip().upper() for item in text.split(",") if item.strip()]
    if not components:
        raise ValueError("--components is required, e.g. A,B,C")
    return components


def load_reference_table(path: Path) -> tuple[dict[str, str], dict[str, float]]:
    """Read columns element,stable_ref_func,orb_reference_eV_atom."""
    stable_ref_funcs: dict[str, str] = {}
    orb_refs: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            el = row["element"].strip().upper()
            stable_ref_funcs[el] = row["stable_ref_func"].strip()
            orb_refs[el] = float(row["orb_reference_eV_atom"])
    return stable_ref_funcs, orb_refs


def present_elements(counts: dict[str, float], components: list[str]) -> list[str]:
    return [el for el in components if counts.get(el, 0.0) > 0]


def site_ratios(counts: dict[str, float], components: list[str]) -> str:
    return " ".join(f"{counts[el]:g}" for el in present_elements(counts, components))


def constituent_line(counts: dict[str, float], components: list[str]) -> str:
    return ":".join(present_elements(counts, components))


def formation_energy_j_mol_atom(
    counts: dict[str, float],
    energy_eV_atom: float,
    components: list[str],
    orb_refs: dict[str, float],
) -> float:
    total = sum(counts.values())
    missing = [el for el in present_elements(counts, components) if el not in orb_refs]
    if missing:
        raise ValueError(f"Missing ORB references for: {', '.join(missing)}")
    ref = sum(counts.get(el, 0.0) / total * orb_refs[el] for el in components)
    return (energy_eV_atom - ref) * EV_PER_ATOM_TO_J_PER_MOL_ATOM


def config_entropy_j_mol_atom(source_cif: str) -> float:
    """Return ideal site-mixing entropy from CIF occupancies, per mol-atom."""
    path = Path(source_cif)
    if not path.exists():
        return 0.0
    structure = Structure.from_file(str(path))
    site_entropy_sum = 0.0
    for site in structure:
        occupancies = [float(v) for v in site.species.as_dict().values()]
        total_occ = sum(occupancies)
        if total_occ <= 0:
            continue
        probabilities = [occ / total_occ for occ in occupancies if occ > 1e-12]
        site_entropy_sum += -sum(p * math.log(p) for p in probabilities)
    if site_entropy_sum <= 1e-14:
        return 0.0
    return R_J_MOL_K * site_entropy_sum / float(structure.composition.num_atoms)


def formula_g_expression(
    counts: dict[str, float],
    formation_j_mol_atom: float,
    config_entropy_j_mol_atom_value: float,
    components: list[str],
    stable_ref_funcs: dict[str, str],
) -> str:
    total = sum(counts.values())
    pieces: list[str] = []
    for el in present_elements(counts, components):
        if el not in stable_ref_funcs:
            raise ValueError(f"Missing stable reference function for {el}")
        pieces.append(f"{counts[el]:g}*{stable_ref_funcs[el]}")
    pieces.append(f"{formation_j_mol_atom * total:.8f}")
    if config_entropy_j_mol_atom_value > 0.0:
        pieces.append(f"-T*{config_entropy_j_mol_atom_value * total:.8f}")
    return " + ".join(pieces)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-tdb", required=True)
    parser.add_argument("--components", required=True, help="Comma-separated components, e.g. A,B,C")
    parser.add_argument(
        "--references",
        required=True,
        help="CSV with columns element,stable_ref_func,orb_reference_eV_atom",
    )
    parser.add_argument("--energies", nargs="+", default=["compound_orb/compound_energies.csv"])
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", default="compound_orb/compound_formation_energies.csv")
    parser.add_argument("--temperature", type=float, required=True)
    args = parser.parse_args()

    components = parse_components(args.components)
    stable_ref_funcs, orb_refs = load_reference_table(Path(args.references))
    base = Path(args.base_tdb).read_text(encoding="utf-8")
    additions: list[str] = [
        "",
        "$ Stoichiometric compound phases from MLIP-relaxed CIF structures.",
        "$ Compound G = stable-element references + MLIP formation energy.",
        "$ Formation energies are 0 K MLIP values with ideal configurational entropy from CIF occupancies.",
    ]
    summary_rows: list[dict[str, object]] = []

    best_rows: dict[str, dict[str, object]] = {}
    for energy_path in args.energies:
        with Path(energy_path).open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                counts = {k.upper(): float(v) for k, v in ast.literal_eval(row["formula_counts"]).items()}
                phase = row["phase"].upper()
                energy_eV_atom = float(row["energy_eV_atom"])
                formation = formation_energy_j_mol_atom(counts, energy_eV_atom, components, orb_refs)
                entropy = config_entropy_j_mol_atom(row["source_cif"])
                free_energy = formation - args.temperature * entropy
                candidate = {
                    "row": row,
                    "counts": counts,
                    "phase": phase,
                    "formation": formation,
                    "config_entropy": entropy,
                    "free_energy": free_energy,
                    "energy_path": energy_path,
                }
                if phase not in best_rows or free_energy < float(best_rows[phase]["free_energy"]):
                    best_rows[phase] = candidate

    for candidate in sorted(best_rows.values(), key=lambda item: str(item["phase"])):
        row = candidate["row"]
        counts = {k.upper(): float(v) for k, v in ast.literal_eval(row["formula_counts"]).items()}
        phase = str(candidate["phase"])
        formation = float(candidate["formation"])
        config_entropy = float(candidate["config_entropy"])
        free_energy = float(candidate["free_energy"])
        energy_eV_atom = float(row["energy_eV_atom"])
        constituents = constituent_line(counts, components)
        g_expr = formula_g_expression(counts, formation, config_entropy, components, stable_ref_funcs)
        additions.extend(
            [
                f" PHASE {phase} % {len(present_elements(counts, components))} {site_ratios(counts, components)} !",
                f"    CONSTITUENT {phase} :{constituents}:!",
                f"   PARAMETER G({phase},{constituents};0) 298.15 {g_expr}; 10000 N REFDUM !",
            ]
        )
        summary_rows.append(
            {
                "phase": phase,
                "display_formula": row["display_formula"],
                "source_cif": row["source_cif"],
                "source_table": candidate["energy_path"],
                "formula_counts": json.dumps(counts, sort_keys=True),
                "energy_eV_atom": energy_eV_atom,
                "formation_eV_atom": formation / EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                "formation_J_mol_atom": formation,
                "config_entropy_J_mol_atom_K": config_entropy,
                f"gibbs_{args.temperature:g}K_eV_atom": free_energy / EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                f"gibbs_{args.temperature:g}K_J_mol_atom": free_energy,
            }
        )

    Path(args.output).write_text(base.rstrip() + "\n" + "\n".join(additions) + "\n", encoding="utf-8")
    with Path(args.summary).open("w", newline="", encoding="utf-8") as f:
        fieldnames = sorted({key for row in summary_rows for key in row})
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.summary}")


if __name__ == "__main__":
    main()
