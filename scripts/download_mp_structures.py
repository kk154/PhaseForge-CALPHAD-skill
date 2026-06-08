#!/usr/bin/env python3
"""Download Materials Project structures for candidate phase compounds."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path


def require_mp_api_key() -> str:
    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        raise SystemExit("MP_API_KEY is not set. Put the key in the environment, not in this script.")
    return api_key


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def get_field(row: dict[str, str], *names: str) -> str:
    lowered = {key.lower(): value for key, value in row.items()}
    for name in names:
        value = lowered.get(name.lower(), "")
        if value:
            return value.strip()
    return ""


def material_id_value(doc: object) -> str:
    value = getattr(doc, "material_id", "")
    return str(value)


def formula_value(doc: object) -> str:
    for attr in ("formula_pretty", "formula_pretty_anonymous", "formula"):
        value = getattr(doc, attr, "")
        if value:
            return str(value)
    return ""


def hull_value(doc: object) -> float | str:
    value = getattr(doc, "energy_above_hull", "")
    try:
        return float(value)
    except Exception:
        return ""


def choose_doc_by_formula(mpr, formula: str, e_above_hull_max: float | None):
    kwargs = {
        "formula": formula,
        "fields": ["material_id", "formula_pretty", "energy_above_hull", "structure"],
    }
    if e_above_hull_max is not None:
        kwargs["energy_above_hull"] = (0, e_above_hull_max)
    docs = list(mpr.materials.summary.search(**kwargs))
    if not docs:
        return None
    return sorted(docs, key=lambda doc: (float(getattr(doc, "energy_above_hull", 9999) or 9999), material_id_value(doc)))[0]


def choose_doc_by_material_id(mpr, material_id: str):
    docs = list(
        mpr.materials.summary.search(
            material_ids=[material_id],
            fields=["material_id", "formula_pretty", "energy_above_hull", "structure"],
        )
    )
    return docs[0] if docs else None


def write_cif(structure, path: Path) -> None:
    from pymatgen.io.cif import CifWriter

    writer = CifWriter(structure)
    writer.write_file(str(path))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidates",
        required=True,
        help="CSV with formula column and optional phase,material_id,source,notes columns",
    )
    parser.add_argument("--output-dir", default="mp_structures")
    parser.add_argument("--manifest", default="mp_structure_manifest.csv")
    parser.add_argument("--e-above-hull-max", type=float, default=0.10)
    parser.add_argument("--allow-metastable-any", action="store_true")
    args = parser.parse_args()

    try:
        from mp_api.client import MPRester
    except Exception as exc:
        raise SystemExit(
            "mp-api is not installed in this Python environment. Run in an environment with mp-api and pymatgen."
        ) from exc

    api_key = require_mp_api_key()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    emax = None if args.allow_metastable_any else args.e_above_hull_max

    with Path(args.candidates).open(newline="", encoding="utf-8") as f:
        candidate_rows = list(csv.DictReader(f))

    manifest_rows: list[dict[str, object]] = []
    with MPRester(api_key) as mpr:
        for row in candidate_rows:
            formula = get_field(row, "formula", "composition", "compound")
            phase = get_field(row, "phase", "phase_name") or formula
            requested_mid = get_field(row, "material_id", "mp_id", "mp-id")
            if not formula and not requested_mid:
                print(f"Skipping row without formula/material_id: {row}", file=sys.stderr)
                continue

            if requested_mid:
                doc = choose_doc_by_material_id(mpr, requested_mid)
            else:
                doc = choose_doc_by_formula(mpr, formula, emax)

            if doc is None:
                manifest_rows.append(
                    {
                        "phase": phase,
                        "formula": formula,
                        "requested_material_id": requested_mid,
                        "material_id": "",
                        "mp_formula": "",
                        "energy_above_hull": "",
                        "cif_path": "",
                        "status": "not_found",
                        "source": get_field(row, "source"),
                        "notes": get_field(row, "notes"),
                    }
                )
                continue

            mid = material_id_value(doc)
            mp_formula = formula_value(doc)
            cif_path = output_dir / f"{safe_name(phase or mp_formula)}_{mid}_{safe_name(mp_formula)}.cif"
            write_cif(getattr(doc, "structure"), cif_path)
            manifest_rows.append(
                {
                    "phase": phase,
                    "formula": formula,
                    "requested_material_id": requested_mid,
                    "material_id": mid,
                    "mp_formula": mp_formula,
                    "energy_above_hull": hull_value(doc),
                    "cif_path": str(cif_path),
                    "status": "downloaded",
                    "source": get_field(row, "source"),
                    "notes": get_field(row, "notes"),
                }
            )
            print(f"Downloaded {phase or formula}: {mid} {mp_formula} -> {cif_path}")

    with Path(args.manifest).open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "phase",
            "formula",
            "requested_material_id",
            "material_id",
            "mp_formula",
            "energy_above_hull",
            "cif_path",
            "status",
            "source",
            "notes",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {args.manifest}")


if __name__ == "__main__":
    main()
