#!/usr/bin/env python3
"""Sample a ternary isothermal section from a TDB with pycalphad."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v


def phase_name(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    return str(value).strip()


def parse_components(text: str) -> list[str]:
    comps = [item.strip().upper() for item in text.split(",") if item.strip()]
    if len(comps) != 3:
        raise ValueError("--components must contain exactly three real components")
    return comps


def active_phases_at(
    dbf: Database,
    components: list[str],
    phases: list[str],
    temperature: float,
    pressure: float,
    xa: float,
    xb: float,
    tol: float,
) -> tuple[str, str]:
    calc_components = components + ["VA"]
    conds = {v.T: temperature, v.P: pressure, v.X(components[0]): xa, v.X(components[1]): xb}
    eq = equilibrium(dbf, calc_components, phases, conds, output="GM")
    active: list[tuple[str, float]] = []
    for phase, amount in zip(eq.Phase.values.ravel(), eq.NP.values.ravel()):
        name = phase_name(phase)
        if name and name != "nan" and np.isfinite(amount) and amount > tol:
            active.append((name, float(amount)))
    if not active:
        return "UNRESOLVED", "UNRESOLVED"
    active.sort(key=lambda item: item[0])
    assemblage = "+".join(sorted({name for name, _ in active}))
    dominant = max(active, key=lambda item: item[1])[0]
    return assemblage, dominant


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tdb", required=True)
    parser.add_argument("--components", required=True, help="Three real components, e.g. ER,FE,TI")
    parser.add_argument("--temperature", type=float, required=True, help="Temperature in K")
    parser.add_argument("--pressure", type=float, default=101325.0)
    parser.add_argument("--step", type=float, default=0.005)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--mode", choices=["assemblage", "dominant"], default="assemblage")
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--exclude-phases", default="", help="Comma-separated phase names to skip")
    args = parser.parse_args()

    components = parse_components(args.components)
    n = int(round(1.0 / args.step))
    if abs(1.0 / n - args.step) > 1e-8:
        raise ValueError("--step must evenly divide 1.0, e.g. 0.05, 0.025, 0.02, 0.005")

    dbf = Database(args.tdb)
    excluded = {item.strip().upper() for item in args.exclude_phases.split(",") if item.strip()}
    phases = sorted(name for name in dbf.phases.keys() if name.upper() not in excluded)

    rows: list[dict[str, object]] = []
    labels: list[str] = []
    for i in range(n + 1):
        xa = i / n
        for j in range(n + 1 - i):
            xb = j / n
            xc = 1.0 - xa - xb
            assemblage, dominant = active_phases_at(
                dbf, components, phases, args.temperature, args.pressure, xa, xb, args.tol
            )
            plot_label = assemblage if args.mode == "assemblage" else dominant
            labels.append(plot_label)
            rows.append(
                {
                    "i": i,
                    "j": j,
                    f"x_{components[0].lower()}": xa,
                    f"x_{components[1].lower()}": xb,
                    f"x_{components[2].lower()}": xc,
                    "component_a": components[0],
                    "component_b": components[1],
                    "component_c": components[2],
                    "temperature_K": args.temperature,
                    "assemblage": assemblage,
                    "dominant": dominant,
                    "plot_label": plot_label,
                }
            )

    with Path(args.csv).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.csv}")
    print("Top labels:")
    for label, count in Counter(labels).most_common(20):
        print(f"  {label}: {count}")


if __name__ == "__main__":
    main()
