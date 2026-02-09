#!/usr/bin/env python3
"""Plot Coulombic energy versus distance from a two-column data file."""

import argparse
import math
from pathlib import Path
from typing import List, Tuple
import matplotlib.pyplot as plt
import prettypyplot as pplt  

def apply_pretty_style() -> None:
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["savefig.edgecolor"] = "white"


def load_distance_energy(path: Path) -> Tuple[List[float], List[float], int]:
    """Load finite distance/energy pairs, skipping comments and bad rows."""
    distances: List[float] = []
    energies: List[float] = []
    dropped = 0

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                dropped += 1
                continue

            try:
                d = float(parts[0])
                e = float(parts[1])
            except ValueError:
                dropped += 1
                continue

            if not (math.isfinite(d) and math.isfinite(e)):
                dropped += 1
                continue

            distances.append(d)
            energies.append(e)

    return distances, energies, dropped


def percentile(sorted_vals: List[float], q: float) -> float:
    """Return percentile from a sorted list using linear interpolation."""
    if not sorted_vals:
        raise ValueError("percentile() received an empty list")

    if q <= 0:
        return sorted_vals[0]
    if q >= 1:
        return sorted_vals[-1]

    pos = (len(sorted_vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    w = pos - lo
    return sorted_vals[lo] * (1.0 - w) + sorted_vals[hi] * w


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot distance_vs_e_energy.dat as Energy vs Distance."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="distance_vs_e_energy.dat",
        help="Input data file (default: distance_vs_e_energy.dat)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="distance_vs_e_energy.png",
        help="Output figure path (default: distance_vs_e_energy.png)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    distances, energies, dropped = load_distance_energy(input_path)
    if not distances:
        raise ValueError(
            f"No finite data points found in {input_path}. "
            "Check for NaN/Inf values in the file."
        )

    apply_pretty_style()

    pairs = sorted(zip(distances, energies), key=lambda t: t[0])
    x_cap = 30.0
    pairs = [p for p in pairs if p[0] <= x_cap]
    if not pairs:
        raise ValueError(f"No points with distance <= {x_cap}")
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    ys_sorted = sorted(ys)

    # Robust visualization window: ignore extreme tails but keep requested top cap at 10.
    y_low = percentile(ys_sorted, 0.01)
    y_high = min(10.0, percentile(ys_sorted, 0.99))
    if y_low >= y_high:
        y_low, y_high = min(ys), 10.0

    inlier_x: List[float] = []
    inlier_y: List[float] = []
    outlier_count = 0
    for x, y in pairs:
        if y_low <= y <= y_high:
            inlier_x.append(x)
            inlier_y.append(y)
        else:
            outlier_count += 1

    fig, ax = plt.subplots(figsize=(9, 5.6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.scatter(xs, ys, s=10, alpha=0.08, color="#4a4a4a", label="All points")
    if inlier_x:
        ax.plot(inlier_x, inlier_y, color="#1f77b4", linewidth=1.8, label="In-range trend")
        ax.scatter(inlier_x, inlier_y, s=12, alpha=0.7, color="#1f77b4")
    ax.set_title("Distance vs Coulombic Electrostatic Energy")
    ax.set_xlabel("Distance")
    ax.set_ylabel("Coulombic Electrostatic Energy")
    ax.set_ylim(y_low, y_high)
    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax.set_xlim(left=0.0, right=x_cap)
    ax.grid(True, alpha=0.25, linestyle="--", linewidth=0.8)
    ax.legend(loc="best", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, facecolor="white", edgecolor="white")

    print(f"Saved figure: {output_path}")
    print(f"Loaded points: {len(distances)}")
    if dropped:
        print(f"Dropped rows (invalid or non-finite): {dropped}")


if __name__ == "__main__":
    main()
