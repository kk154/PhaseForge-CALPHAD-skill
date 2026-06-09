#!/usr/bin/env python3
"""Fit a practical CALPHAD LIQUID model from SQS/LAMMPS liquid enthalpies."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


EV_PER_ATOM_TO_J_PER_MOL_ATOM = 96485.33212331002

MELTING_POINTS_K = {
    "H": 14.01,
    "HE": 0.0,
    "LI": 453.65,
    "BE": 1560.0,
    "B": 2349.0,
    "C": 3800.0,
    "N": 63.15,
    "O": 54.36,
    "F": 53.53,
    "NE": 24.56,
    "NA": 370.87,
    "MG": 923.0,
    "AL": 933.47,
    "SI": 1687.0,
    "P": 317.3,
    "S": 388.36,
    "CL": 171.6,
    "AR": 83.8,
    "K": 336.53,
    "CA": 1115.0,
    "SC": 1814.0,
    "TI": 1941.0,
    "V": 2183.0,
    "CR": 2180.0,
    "MN": 1519.0,
    "FE": 1811.0,
    "CO": 1768.0,
    "NI": 1728.0,
    "CU": 1357.77,
    "ZN": 692.88,
    "GA": 302.91,
    "GE": 1211.4,
    "AS": 1090.0,
    "SE": 453.0,
    "BR": 265.8,
    "KR": 115.79,
    "RB": 312.46,
    "SR": 1050.0,
    "Y": 1799.0,
    "ZR": 2128.0,
    "NB": 2750.0,
    "MO": 2896.0,
    "TC": 2430.0,
    "RU": 2607.0,
    "RH": 2237.0,
    "PD": 1828.05,
    "AG": 1234.93,
    "CD": 594.22,
    "IN": 429.75,
    "SN": 505.08,
    "SB": 903.78,
    "I": 386.85,
    "TE": 722.66,
    "XE": 161.4,
    "CS": 301.59,
    "BA": 1000.0,
    "LA": 1193.0,
    "CE": 1068.0,
    "PR": 1208.0,
    "ND": 1297.0,
    "PM": 1315.0,
    "SM": 1345.0,
    "EU": 1099.0,
    "GD": 1585.0,
    "TB": 1629.0,
    "DY": 1680.0,
    "HO": 1734.0,
    "ER": 1802.0,
    "TM": 1818.0,
    "YB": 1097.0,
    "LU": 1925.0,
    "HF": 2506.0,
    "TA": 3290.0,
    "W": 3695.0,
    "RE": 3459.0,
    "OS": 3306.0,
    "IR": 2719.0,
    "PT": 2041.4,
    "AU": 1337.33,
    "HG": 234.43,
    "TL": 577.0,
    "PB": 600.61,
    "BI": 544.7,
    "PO": 527.0,
    "AT": 575.0,
    "RN": 202.0,
    "FR": 300.0,
    "RA": 973.0,
    "AC": 1323.0,
    "TH": 2115.0,
    "PA": 1841.0,
    "U": 1405.3,
    "NP": 917.0,
    "PU": 912.5,
    "AM": 1449.0,
    "CM": 1613.0,
    "BK": 1259.0,
    "CF": 1173.0,
    "ES": 1133.0,
    "FM": 1125.0,
    "MD": 1100.0,
    "NO": 1100.0,
    "LR": 1900.0,
    "RF": 2400.0,
}


def parse_components(text: str) -> list[str]:
    components = [item.strip().upper() for item in text.split(",") if item.strip()]
    if len(components) != 3:
        raise ValueError("--components must contain exactly three elements")
    return components


def load_reference_table(path: Path) -> tuple[dict[str, str], dict[str, float]]:
    stable_ref_funcs: dict[str, str] = {}
    solid_refs: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            el = row["element"].strip().upper()
            stable_ref_funcs[el] = row["stable_ref_func"].strip()
            solid_refs[el] = float(row["orb_reference_eV_atom"])
    return stable_ref_funcs, solid_refs


def parse_composition(row: dict[str, str]) -> dict[str, float]:
    if row.get("composition_json"):
        return {str(k).upper(): float(v) for k, v in json.loads(row["composition_json"]).items()}
    comp: dict[str, float] = {}
    for part in row["composition"].split(";"):
        if not part:
            continue
        el, value = part.split(":", 1)
        comp[el.strip().upper()] = float(value)
    total = sum(comp.values())
    return {el: value / total for el, value in comp.items()}


def is_pure(comp: dict[str, float], element: str) -> bool:
    return comp.get(element, 0.0) > 1.0 - 1e-8


def composition_key(comp: dict[str, float], components: list[str]) -> tuple[float, ...]:
    return tuple(round(comp.get(el, 0.0), 10) for el in components)


def load_liquid_rows(path: Path, components: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            comp = parse_composition(row)
            if any(el not in components for el in comp):
                raise ValueError(f"Liquid composition contains element outside --components: {comp}")
            rows.append(
                {
                    "composition": comp,
                    "key": composition_key(comp, components),
                    "temperature_K": float(row["temperature_K"]),
                    "enthalpy_eV_atom": float(row["enthalpy_eV_atom"]),
                    "source_dir": row.get("source_dir", ""),
                }
            )
    return rows


def fit_linear_t(values: list[tuple[float, float]]) -> tuple[float, float]:
    if not values:
        raise ValueError("No values to fit")
    if len(values) == 1:
        return values[0][1], 0.0
    temps = np.asarray([item[0] for item in values], dtype=float)
    ys = np.asarray([item[1] for item in values], dtype=float)
    b, a = np.polyfit(temps, ys, 1)
    return float(a), float(b)


def parameter_expr(a_j: float, b_j_per_k: float) -> str:
    if abs(b_j_per_k) < 1e-10:
        return f"{a_j:.8f}"
    sign = "+" if b_j_per_k >= 0 else "-"
    return f"{a_j:.8f} {sign} {abs(b_j_per_k):.8f}*T"


def strip_liquid_from_tdb(text: str) -> str:
    lines = text.splitlines()
    output: list[str] = []
    skipping = False
    for line in lines:
        upper = line.upper()
        if re.match(r"\s*PHASE\s+LIQUID\b", upper):
            continue
        if re.match(r"\s*CONSTITUENT\s+LIQUID\b", upper):
            continue
        if re.match(r"\s*PARAMETER\s+[GL]\(LIQUID,", upper):
            continue
        output.append(line)
    return "\n".join(output).rstrip() + "\n"


def pure_liquid_enthalpies(rows: list[dict[str, object]], components: list[str]) -> dict[str, list[tuple[float, float]]]:
    result: dict[str, list[tuple[float, float]]] = {el: [] for el in components}
    for row in rows:
        comp = row["composition"]
        assert isinstance(comp, dict)
        for el in components:
            if is_pure(comp, el):
                result[el].append((float(row["temperature_K"]), float(row["enthalpy_eV_atom"])))
    return result


def interpolated_pure_enthalpy(pure_fits: dict[str, tuple[float, float]], element: str, temperature: float) -> float:
    a, b = pure_fits[element]
    return a + b * temperature


def excess_enthalpy(row: dict[str, object], components: list[str], pure_fits: dict[str, tuple[float, float]]) -> float:
    comp = row["composition"]
    assert isinstance(comp, dict)
    temperature = float(row["temperature_K"])
    reference = sum(comp.get(el, 0.0) * interpolated_pure_enthalpy(pure_fits, el, temperature) for el in components)
    return float(row["enthalpy_eV_atom"]) - reference


def fit_binary_excess(
    rows: list[dict[str, object]],
    components: list[str],
    pure_fits: dict[str, tuple[float, float]],
) -> list[dict[str, object]]:
    parameters: list[dict[str, object]] = []
    for i, a_el in enumerate(components):
        for b_el in components[i + 1 :]:
            values: list[tuple[float, float]] = []
            for row in rows:
                comp = row["composition"]
                assert isinstance(comp, dict)
                present = {el for el, x in comp.items() if x > 1e-8}
                if present != {a_el, b_el}:
                    continue
                xa = comp.get(a_el, 0.0)
                xb = comp.get(b_el, 0.0)
                denom = xa * xb
                if denom <= 1e-10:
                    continue
                values.append((float(row["temperature_K"]), excess_enthalpy(row, components, pure_fits) / denom))
            if values:
                l0_eV, l0_slope_eV = fit_linear_t(values)
                parameters.append(
                    {
                        "kind": "binary",
                        "elements": (a_el, b_el),
                        "order": 0,
                        "a_J": l0_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                        "b_J_per_K": l0_slope_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                        "num_points": len(values),
                    }
                )
    return parameters


def fit_ternary_excess(
    rows: list[dict[str, object]],
    components: list[str],
    pure_fits: dict[str, tuple[float, float]],
    binary_params: list[dict[str, object]],
) -> list[dict[str, object]]:
    values: list[tuple[float, float]] = []
    for row in rows:
        comp = row["composition"]
        assert isinstance(comp, dict)
        if any(comp.get(el, 0.0) <= 1e-8 for el in components):
            continue
        temperature = float(row["temperature_K"])
        predicted = 0.0
        for param in binary_params:
            a_el, b_el = param["elements"]
            assert isinstance(a_el, str) and isinstance(b_el, str)
            l0 = (float(param["a_J"]) + float(param["b_J_per_K"]) * temperature) / EV_PER_ATOM_TO_J_PER_MOL_ATOM
            predicted += comp[a_el] * comp[b_el] * l0
        residual = excess_enthalpy(row, components, pure_fits) - predicted
        denom = math.prod(comp[el] for el in components)
        if denom > 1e-12:
            values.append((temperature, residual / denom))
    if not values:
        return []
    a_eV, b_eV = fit_linear_t(values)
    return [
        {
            "kind": "ternary",
            "elements": tuple(components),
            "a_J": a_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
            "b_J_per_K": b_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
            "num_points": len(values),
        }
    ]


def liquid_phase_block(
    components: list[str],
    stable_ref_funcs: dict[str, str],
    solid_refs: dict[str, float],
    pure_fits: dict[str, tuple[float, float]],
    binary_params: list[dict[str, object]],
    ternary_params: list[dict[str, object]],
) -> str:
    lines = [
        "",
        "$ LIQUID phase fitted from SQS/LAMMPS MLIP liquid enthalpies.",
        "$ Pure liquid terms use G_SER + dH_fus*(1-T/Tm); excess terms use L(T)=a+b*T.",
        " PHASE LIQUID % 1 1 !",
        f"    CONSTITUENT LIQUID :{','.join(components)}:!",
    ]
    for el in components:
        if el not in MELTING_POINTS_K or MELTING_POINTS_K[el] <= 0:
            raise ValueError(f"Missing positive melting point for {el}")
        if el not in stable_ref_funcs:
            raise ValueError(f"Missing stable reference function for {el}")
        pure_a, pure_b = pure_fits[el]
        h_liq_at_tm = pure_a + pure_b * MELTING_POINTS_K[el]
        d_h_fus = (h_liq_at_tm - solid_refs[el]) * EV_PER_ATOM_TO_J_PER_MOL_ATOM
        tm = MELTING_POINTS_K[el]
        expr = f"{stable_ref_funcs[el]} + ({d_h_fus:.8f})*(1 - T/{tm:.8f})"
        lines.append(f"   PARAMETER G(LIQUID,{el};0) 298.15 {expr}; 10000 N REFDUM !")
    for param in binary_params:
        a_el, b_el = param["elements"]
        expr = parameter_expr(float(param["a_J"]), float(param["b_J_per_K"]))
        lines.append(f"   PARAMETER G(LIQUID,{a_el},{b_el};0) 298.15 {expr}; 10000 N REFDUM !")
    for param in ternary_params:
        a_el, b_el, c_el = param["elements"]
        expr = parameter_expr(float(param["a_J"]), float(param["b_J_per_K"]))
        lines.append(f"   PARAMETER G(LIQUID,{a_el},{b_el},{c_el};0) 298.15 {expr}; 10000 N REFDUM !")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-tdb", required=True)
    parser.add_argument("--liquid-md", required=True)
    parser.add_argument("--references", required=True, help="CSV with element,stable_ref_func,orb_reference_eV_atom")
    parser.add_argument("--components", required=True)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", default="liquid_fit_summary.csv")
    args = parser.parse_args()

    components = parse_components(args.components)
    stable_ref_funcs, solid_refs = load_reference_table(Path(args.references))
    rows = load_liquid_rows(Path(args.liquid_md), components)
    pure_rows = pure_liquid_enthalpies(rows, components)
    missing_pure = [el for el, values in pure_rows.items() if not values]
    if missing_pure:
        raise ValueError(f"Missing pure-liquid MD rows for: {', '.join(missing_pure)}")

    pure_fits = {el: fit_linear_t(values) for el, values in pure_rows.items()}
    binary_params = fit_binary_excess(rows, components, pure_fits)
    ternary_params = fit_ternary_excess(rows, components, pure_fits, binary_params)
    base = strip_liquid_from_tdb(Path(args.base_tdb).read_text(encoding="utf-8"))
    block = liquid_phase_block(components, stable_ref_funcs, solid_refs, pure_fits, binary_params, ternary_params)
    Path(args.output).write_text(base.rstrip() + "\n" + block, encoding="utf-8")

    with Path(args.summary).open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["kind", "elements", "a_J", "b_J_per_K", "num_points", "target_temperature_K"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for el, (a_eV, b_eV) in pure_fits.items():
            writer.writerow(
                {
                    "kind": "pure",
                    "elements": el,
                    "a_J": a_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                    "b_J_per_K": b_eV * EV_PER_ATOM_TO_J_PER_MOL_ATOM,
                    "num_points": len(pure_rows[el]),
                    "target_temperature_K": args.temperature,
                }
            )
        for param in binary_params + ternary_params:
            writer.writerow(
                {
                    "kind": param["kind"],
                    "elements": ",".join(param["elements"]),
                    "a_J": param["a_J"],
                    "b_J_per_K": param["b_J_per_K"],
                    "num_points": param["num_points"],
                    "target_temperature_K": args.temperature,
                }
            )

    try:
        from pycalphad import Database
    except ImportError:
        print("WARNING: pycalphad is not installed; skipped TDB validation.")
    else:
        try:
            Database(args.output)
        except Exception as exc:  # pragma: no cover - validation depends on optional pycalphad.
            raise SystemExit(f"Wrote {args.output}, but pycalphad validation failed: {exc}") from exc

    print(f"Wrote {args.output}")
    print(f"Wrote {args.summary}")


if __name__ == "__main__":
    main()
