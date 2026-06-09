#!/usr/bin/env python3
"""Run one liquid LAMMPS/MLIP MD calculation and summarize enthalpy."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData


BAR_A3_TO_EV = 6.241509074e-7


def parse_triplet(text: str) -> list[int]:
    values = [int(item) for item in text.split(",")]
    if len(values) != 3:
        raise ValueError("Expected three comma-separated integers")
    return values


def species_labels(structure: Structure, atomlabel_path: Path) -> list[str]:
    if atomlabel_path.exists():
        labels = [item.strip() for item in atomlabel_path.read_text(encoding="utf-8").split() if item.strip()]
        if labels:
            return labels
    labels: list[str] = []
    for site in structure:
        label = str(site.specie)
        if label not in labels:
            labels.append(label)
    return labels


def write_lammps_input(
    path: Path,
    *,
    data_file: str,
    temperature: float,
    timestep: float,
    equil_steps: int,
    prod_steps: int,
    sample_every: int,
    pair_style: str,
    gnnp_root: str,
    mlip: str,
    model: str,
    labels: list[str],
) -> None:
    window = max(1, prod_steps // max(1, sample_every))
    average_every = sample_every * window
    label_text = " ".join(labels)
    path.write_text(
        f"""variable        pe_atom equal pe/count(all)
variable        pv_atom equal press*vol*{BAR_A3_TO_EV:.15g}/count(all)
variable        h_atom equal v_pe_atom+v_pv_atom
variable        p equal press
variable        v equal vol

units           metal
dimension       3
boundary        p p p
atom_style      atomic

read_data       ./{data_file}

pair_style      {pair_style} {gnnp_root}
pair_coeff      * * {mlip} {model} {label_text}

timestep        {timestep}

reset_timestep  0
thermo          100
thermo_style    custom step temp press pe ke etotal vol lx ly lz

min_style       cg
minimize        1e-25 1e-25 5000 10000

reset_timestep  0
velocity        all create {temperature:.10g} 93723 mom yes rot no

fix             EQ all npt temp {temperature:.10g} {temperature:.10g} $(10.0*dt) iso 0.0 0.0 $(100.0*dt)
run             {equil_steps}
unfix           EQ

fix             PROD all npt temp {temperature:.10g} {temperature:.10g} $(10.0*dt) iso 0.0 0.0 $(100.0*dt)
fix             AVG all ave/time {sample_every} {window} {average_every} v_pe_atom v_pv_atom v_h_atom v_p v_v file liquid_averages.dat ave one
fix_modify      AVG title1 "# pe_eV_atom pv_eV_atom enthalpy_eV_atom press_bar volume_A3"
thermo          100
thermo_style    custom step temp press pe ke etotal vol v_pe_atom v_pv_atom v_h_atom
run             {prod_steps}
unfix           PROD
unfix           AVG
""",
        encoding="utf-8",
    )


def parse_average(path: Path) -> dict[str, float]:
    if not path.exists():
        raise FileNotFoundError(f"LAMMPS average file not found: {path}")
    data_line = ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            data_line = stripped
    if not data_line:
        raise ValueError(f"No averaged data rows found in {path}")
    values = [float(item) for item in data_line.split()]
    if len(values) == 6:
        _, pe, pv, h, press, volume = values
    elif len(values) == 5:
        pe, pv, h, press, volume = values
    else:
        raise ValueError(f"Unexpected liquid average row in {path}: {data_line}")
    return {
        "pe_eV_atom": pe,
        "pv_eV_atom": pv,
        "enthalpy_eV_atom": h,
        "press_bar": press,
        "volume_A3": volume,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poscar", default="POSCAR")
    parser.add_argument("--atomlabel", default="../atomlabel.tmp")
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--mlip", default=os.environ.get("PHASEFORGE_MLIP", "ORB"))
    parser.add_argument("--model", default=os.environ.get("PHASEFORGE_MODEL", "orb-v3-conservative-inf-omat"))
    parser.add_argument("--lammps", default=os.environ.get("PHASEFORGE_LMP", "lmp"))
    parser.add_argument("--gnnp-root", default=os.environ.get("PHASEFORGE_GNNP_ROOT", str(Path.home() / ".local/share/ML-GNNP")))
    parser.add_argument("--pair-style", default=os.environ.get("PHASEFORGE_GNNP_PAIR_STYLE", "gnnp/gpu"))
    parser.add_argument("--supercell", default="2,2,2")
    parser.add_argument("--equil-steps", type=int, default=int(os.environ.get("PHASEFORGE_LIQUID_EQUIL_STEPS", "30000")))
    parser.add_argument("--prod-steps", type=int, default=int(os.environ.get("PHASEFORGE_LIQUID_MD_STEPS", "30000")))
    parser.add_argument("--sample-every", type=int, default=10)
    parser.add_argument("--timestep", type=float, default=0.001)
    parser.add_argument("--no-run", action="store_true")
    args = parser.parse_args()

    structure = Structure.from_file(args.poscar)
    labels = species_labels(structure, Path(args.atomlabel))
    supercell = structure * parse_triplet(args.supercell)
    LammpsData.from_structure(supercell, atom_style="atomic").write_file("POSCAR.data")
    write_lammps_input(
        Path("lammps.in"),
        data_file="POSCAR.data",
        temperature=args.temperature,
        timestep=args.timestep,
        equil_steps=args.equil_steps,
        prod_steps=args.prod_steps,
        sample_every=args.sample_every,
        pair_style=args.pair_style,
        gnnp_root=args.gnnp_root,
        mlip=args.mlip,
        model=args.model,
        labels=labels,
    )

    if args.no_run:
        print("LAMMPS execution skipped due to --no-run.")
        return

    with Path("lammps.log").open("w", encoding="utf-8") as log:
        subprocess.run([args.lammps, "-in", "lammps.in"], check=True, stdout=log, stderr=subprocess.STDOUT)

    averages = parse_average(Path("liquid_averages.dat"))
    counts = {str(el).upper(): float(amount) for el, amount in structure.composition.as_dict().items()}
    total = sum(counts.values())
    composition = {el: amount / total for el, amount in counts.items()}
    summary = {
        "phase": "LIQUID",
        "temperature_K": args.temperature,
        "composition": composition,
        "formula_counts": counts,
        "natoms_unitcell": float(len(structure)),
        "supercell": args.supercell,
        "equil_steps": args.equil_steps,
        "prod_steps": args.prod_steps,
        "sample_every": args.sample_every,
        "mlip": args.mlip,
        "model": args.model,
        **averages,
    }
    Path("liquid_md_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    Path("energy").write_text(f"{averages['enthalpy_eV_atom']:.10f}\n", encoding="utf-8")
    print(f"Wrote liquid_md_summary.json and energy={averages['enthalpy_eV_atom']:.10f} eV/atom")


if __name__ == "__main__":
    main()
