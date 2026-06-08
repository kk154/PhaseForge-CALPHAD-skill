#!/usr/bin/env python3
"""Draw a Thermo-Calc-style colored ternary isothermal section."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def barycentric_to_xy(xa: float, xb: float, xc: float) -> tuple[float, float]:
    total = xa + xb + xc
    xa, xb, xc = xa / total, xb / total, xc / total
    return xb + 0.5 * xc, (np.sqrt(3.0) / 2.0) * xc


def split_components(text: str | None) -> list[str] | None:
    if not text:
        return None
    comps = [item.strip().upper() for item in text.split(",") if item.strip()]
    if len(comps) != 3:
        raise ValueError("--components must contain exactly three components")
    return comps


def clean_phase_name(name: str, aliases: dict[str, str]) -> str:
    return aliases.get(name, name)


def clean_assemblage(label: str, aliases: dict[str, str]) -> str:
    parts = [clean_phase_name(part, aliases) for part in label.split("+")]
    text = "+".join(parts)
    if len(text) <= 42:
        return text
    lines: list[str] = []
    current = ""
    for part in parts:
        trial = part if not current else current + "+" + part
        if len(trial) > 36 and current:
            lines.append(current)
            current = part
        else:
            current = trial
    if current:
        lines.append(current)
    return "\n".join(lines)


def read_aliases(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    aliases: dict[str, str] = {}
    with Path(path).open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            aliases[row["phase"]] = row["label"]
    return aliases


def infer_components(fieldnames: list[str], requested: list[str] | None) -> list[str]:
    if requested:
        return requested
    if {"component_a", "component_b", "component_c"}.issubset(fieldnames):
        return []
    x_cols = [name for name in fieldnames if name.startswith("x_")]
    if len(x_cols) < 3:
        raise ValueError("Could not infer composition columns; pass --components A,B,C")
    return [name[2:].upper() for name in x_cols[:3]]


def read_grid(path: Path, label_column: str, components_arg: str | None) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    requested = split_components(components_arg)
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Empty grid: {path}")
        inferred = infer_components(reader.fieldnames, requested)
        for row in reader:
            components = inferred
            if not components:
                components = [row["component_a"].upper(), row["component_b"].upper(), row["component_c"].upper()]
            xa = float(row[f"x_{components[0].lower()}"])
            xb = float(row[f"x_{components[1].lower()}"])
            xc = float(row[f"x_{components[2].lower()}"])
            x, y = barycentric_to_xy(xa, xb, xc)
            rows.append(
                {
                    "i": int(row["i"]) if row.get("i", "") != "" else -1,
                    "j": int(row["j"]) if row.get("j", "") != "" else -1,
                    "xa": xa,
                    "xb": xb,
                    "xc": xc,
                    "x": x,
                    "y": y,
                    "label": row[label_column],
                    "dominant": row.get("dominant", row[label_column]),
                    "assemblage": row.get("assemblage", row[label_column]),
                }
            )
    if not rows:
        raise ValueError(f"No rows read from {path}")
    return rows, components


def draw_ternary_grid(ax: plt.Axes, step: float = 0.1) -> None:
    vals = np.arange(step, 1.0, step)
    for value in vals:
        grid_lines = [
            (barycentric_to_xy(value, 0, 1 - value), barycentric_to_xy(value, 1 - value, 0)),
            (barycentric_to_xy(0, value, 1 - value), barycentric_to_xy(1 - value, value, 0)),
            (barycentric_to_xy(0, 1 - value, value), barycentric_to_xy(1 - value, 0, value)),
        ]
        for p1, p2 in grid_lines:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.82", linewidth=0.45, zorder=1)


def draw_ticks(ax: plt.Axes, components: list[str], step: float = 0.1) -> None:
    vals = np.arange(0, 1.01, step)
    for value in vals:
        x, y = barycentric_to_xy(1 - value, value, 0)
        ax.text(x, y - 0.026, f"{value:.1f}", ha="center", va="top", fontsize=7)
        x, y = barycentric_to_xy(1 - value, 0, value)
        ax.text(x - 0.024, y, f"{value:.1f}", ha="right", va="center", fontsize=7, rotation=60)
        x, y = barycentric_to_xy(0, 1 - value, value)
        ax.text(x + 0.024, y, f"{value:.1f}", ha="left", va="center", fontsize=7, rotation=-60)
    ax.text(0.5, -0.066, f"Mole fraction {components[1]}", ha="center", va="top", fontsize=9)


def boundary_segments(rows: list[dict[str, object]]) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    xs = np.array([float(row["x"]) for row in rows])
    ys = np.array([float(row["y"]) for row in rows])
    labels = [str(row["label"]) for row in rows]
    tri = mtri.Triangulation(xs, ys)
    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for triangle in tri.triangles:
        tri_labels = [labels[int(idx)] for idx in triangle]
        if len(set(tri_labels)) == 1:
            continue
        points = [(xs[int(idx)], ys[int(idx)]) for idx in triangle]
        mids: list[tuple[float, float]] = []
        for a, b in [(0, 1), (1, 2), (2, 0)]:
            if tri_labels[a] != tri_labels[b]:
                mids.append((0.5 * (points[a][0] + points[b][0]), 0.5 * (points[a][1] + points[b][1])))
        if len(set(tri_labels)) == 2 and len(mids) == 2:
            segments.append((mids[0], mids[1]))
        elif len(mids) >= 3:
            center = (float(np.mean([p[0] for p in points])), float(np.mean([p[1] for p in points])))
            segments.extend((center, mid) for mid in mids)
    return segments


def straight_boundary_segments(
    rows: list[dict[str, object]],
    min_points: int,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    by_index = {
        (int(row["i"]), int(row["j"])): row
        for row in rows
        if int(row.get("i", -1)) >= 0 and int(row.get("j", -1)) >= 0
    }
    groups: dict[tuple[str, str], list[tuple[float, float]]] = {}
    for (i, j), row in by_index.items():
        label = str(row["label"])
        if label == "UNRESOLVED":
            continue
        for di, dj in [(1, 0), (0, 1), (1, -1)]:
            other = by_index.get((i + di, j + dj))
            if other is None:
                continue
            other_label = str(other["label"])
            if other_label == "UNRESOLVED" or other_label == label:
                continue
            pair = tuple(sorted((label, other_label)))
            groups.setdefault(pair, []).append(
                (
                    0.5 * (float(row["x"]) + float(other["x"])),
                    0.5 * (float(row["y"]) + float(other["y"])),
                )
            )

    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for points in groups.values():
        if len(points) < min_points:
            continue
        data = np.asarray(points, dtype=float)
        center = np.mean(data, axis=0)
        centered = data - center
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        direction = vt[0]
        projections = centered @ direction
        p1 = center + direction * float(np.min(projections))
        p2 = center + direction * float(np.max(projections))
        segments.append(((float(p1[0]), float(p1[1])), (float(p2[0]), float(p2[1]))))
    return segments


def triangle_face_labels(rows: list[dict[str, object]], tri: mtri.Triangulation) -> list[str]:
    labels = [str(row["label"]) for row in rows]
    face_labels: list[str] = []
    for triangle in tri.triangles:
        values = [labels[int(idx)] for idx in triangle]
        face_labels.append(Counter(values).most_common(1)[0][0])
    return face_labels


def label_centroids(rows: list[dict[str, object]], max_labels: int) -> list[tuple[str, float, float, int]]:
    grouped: dict[str, list[tuple[float, float]]] = {}
    for row in rows:
        label = str(row["label"])
        if label == "UNRESOLVED":
            continue
        grouped.setdefault(label, []).append((float(row["x"]), float(row["y"])))
    centroids = []
    for label, points in grouped.items():
        if len(points) < 2:
            continue
        centroids.append(
            (
                label,
                float(np.median([p[0] for p in points])),
                float(np.median([p[1] for p in points])),
                len(points),
            )
        )
    centroids.sort(key=lambda item: item[3], reverse=True)
    return centroids[:max_labels]


def read_markers(path: str | None) -> list[tuple[str, tuple[float, float, float]]]:
    if not path:
        return []
    markers = []
    with Path(path).open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            markers.append((row["label"], (float(row["x_a"]), float(row["x_b"]), float(row["x_c"]))))
    return markers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--components", help="Three components in vertex order A,B,C; inferred from grid if omitted")
    parser.add_argument("--temperature", help="Temperature label, e.g. 1173.15 K or 900 C")
    parser.add_argument("--title")
    parser.add_argument("--label-column", choices=["assemblage", "dominant", "plot_label"], default="assemblage")
    parser.add_argument("--phase-aliases", help="CSV with columns phase,label")
    parser.add_argument("--markers", help="CSV with columns label,x_a,x_b,x_c")
    parser.add_argument("--max-internal-labels", type=int, default=18)
    parser.add_argument("--straight-boundaries", action="store_true")
    parser.add_argument("--min-straight-boundary-points", type=int, default=25)
    args = parser.parse_args()

    aliases = read_aliases(args.phase_aliases)
    rows, components = read_grid(Path(args.grid), args.label_column, args.components)
    xs = np.array([float(row["x"]) for row in rows])
    ys = np.array([float(row["y"]) for row in rows])
    tri = mtri.Triangulation(xs, ys)
    face_labels = triangle_face_labels(rows, tri)
    unique = sorted(set(face_labels), key=lambda label: (-face_labels.count(label), label))

    base_colors = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    colors = {label: base_colors[idx % len(base_colors)] for idx, label in enumerate(unique)}
    label_to_index = {label: idx for idx, label in enumerate(unique)}
    face_values = np.array([label_to_index[label] for label in face_labels], dtype=float)
    cmap = ListedColormap([colors[label] for label in unique])

    fig, ax = plt.subplots(figsize=(13.2, 8.2), constrained_layout=True)
    draw_ternary_grid(ax)
    ax.tripcolor(tri, facecolors=face_values, cmap=cmap, edgecolors="none", alpha=0.62, zorder=0)

    if args.straight_boundaries:
        segments = straight_boundary_segments(rows, args.min_straight_boundary_points)
    else:
        segments = boundary_segments(rows)
    for (x1, y1), (x2, y2) in segments:
        ax.plot([x1, x2], [y1, y2], color="black", linewidth=0.65, zorder=3)

    triangle = np.array(
        [
            barycentric_to_xy(1, 0, 0),
            barycentric_to_xy(0, 1, 0),
            barycentric_to_xy(0, 0, 1),
            barycentric_to_xy(1, 0, 0),
        ]
    )
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1.6, zorder=4)
    draw_ticks(ax, components)

    for label, x, y, _count in label_centroids(rows, args.max_internal_labels):
        text = clean_assemblage(label, aliases)
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=6.5,
            color="black",
            bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 1.2},
            zorder=5,
        )

    for name, comp in read_markers(args.markers):
        x, y = barycentric_to_xy(*comp)
        ax.plot([x], [y], marker="o", color="black", markersize=4.0, zorder=6)
        ax.text(x + 0.012, y + 0.010, name, fontsize=7.5, zorder=6)

    ax.text(-0.045, -0.045, components[0], ha="right", va="top", fontsize=15)
    ax.text(1.045, -0.045, components[1], ha="left", va="top", fontsize=15)
    ax.text(0.5, np.sqrt(3) / 2 + 0.026, components[2], ha="center", va="bottom", fontsize=15)
    if args.title:
        title = args.title
    else:
        temp = f" at {args.temperature}" if args.temperature else ""
        title = f"{'-'.join(components)} isothermal section{temp} - Thermo-Calc style"
    ax.set_title(title, fontsize=14)

    legend_handles = [
        Patch(facecolor=colors[label], edgecolor="black", linewidth=0.25, label=clean_assemblage(label, aliases))
        for label in unique
    ]
    ax.legend(
        handles=legend_handles,
        title="Phase assemblage",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=5.8,
        title_fontsize=8,
        frameon=True,
        borderpad=0.7,
        labelspacing=0.35,
    )

    ax.set_xlim(-0.075, 1.075)
    ax.set_ylim(-0.09, np.sqrt(3) / 2 + 0.09)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(args.output, dpi=450)
    pdf_path = str(Path(args.output).with_suffix(".pdf"))
    fig.savefig(pdf_path)
    print(f"Wrote {args.output}")
    print(f"Wrote {pdf_path}")
    print(f"Assemblages: {len(unique)}")
    if args.straight_boundaries:
        print(f"Straight boundary segments: {len(segments)}")


if __name__ == "__main__":
    main()
