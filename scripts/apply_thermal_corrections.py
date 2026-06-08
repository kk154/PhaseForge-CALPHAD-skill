#!/usr/bin/env python3
"""Apply finite-temperature Gibbs corrections to selected TDB phases.

The correction CSV must contain:
phase,correction_eV_atom

or:
phase,correction_J_mol_atom

Corrections are applied as constant Gibbs shifts per mol-atom. For
multi-sublattice phases, each existing G endmember parameter is shifted by
correction * sum(site ratios), because pycalphad reports GM per mol-atom after
dividing by the site-ratio sum.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002


def parse_phase_site_totals(tdb_text: str) -> dict[str, float]:
    totals: dict[str, float] = {}
    pattern = re.compile(r"^\s*PHASE\s+(\S+)\s+\S+\s+\d+\s+(.+?)\s*!", re.MULTILINE | re.IGNORECASE)
    for match in pattern.finditer(tdb_text):
        phase = match.group(1).upper()
        ratios = [float(token) for token in match.group(2).split()]
        totals[phase] = sum(ratios)
    return totals


def load_corrections(path: Path) -> dict[str, float]:
    corrections: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            phase = row["phase"].upper()
            if row.get("correction_J_mol_atom"):
                value = float(row["correction_J_mol_atom"])
            else:
                value = float(row["correction_eV_atom"]) * EV_PER_ATOM_TO_J_PER_MOL_ATOM
            corrections[phase] = value
    return corrections


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-tdb", required=True)
    parser.add_argument("--corrections", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    text = Path(args.base_tdb).read_text(encoding="utf-8")
    phase_site_totals = parse_phase_site_totals(text)
    corrections = load_corrections(Path(args.corrections))

    for phase in corrections:
        if phase not in phase_site_totals:
            raise ValueError(f"Phase {phase} not found in {args.base_tdb}")

    parameter_pattern = re.compile(
        r"(?P<prefix>\s*PARAMETER\s+G\((?P<phase>[^,]+),(?P<constituents>[^;]+);0\)\s+298\.15\s+)"
        r"(?P<expr>.*?)(?P<suffix>;\s+10000\s+N\s+REFDUM\s*!)",
        re.IGNORECASE,
    )
    applied_counts = {phase: 0 for phase in corrections}

    def replace(match: re.Match[str]) -> str:
        phase = match.group("phase").upper()
        if phase not in corrections:
            return match.group(0)
        shift = corrections[phase] * phase_site_totals[phase]
        applied_counts[phase] += 1
        return f"{match.group('prefix')}({match.group('expr')}) + {shift:.8f}{match.group('suffix')}"

    updated_text = parameter_pattern.sub(replace, text)
    missing = [phase for phase, count in applied_counts.items() if count == 0]
    if missing:
        raise ValueError(f"No G endmember parameters found for: {', '.join(missing)}")

    header = [
        "$ Finite-temperature Gibbs corrections applied to existing G parameters.",
        "$ Values are user-supplied per mol-atom and expanded by phase site-ratio sums.",
    ]
    Path(args.output).write_text("\n".join(header) + "\n" + updated_text, encoding="utf-8")
    print(f"Wrote {args.output}")
    for phase, count in sorted(applied_counts.items()):
        print(f"Applied correction to {count} G parameters for {phase}")


if __name__ == "__main__":
    main()
