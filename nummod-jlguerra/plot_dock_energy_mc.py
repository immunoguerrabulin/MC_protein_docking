#!/usr/bin/env python3
"""Plot MC docking energy trace from dock_energy_mc.dat."""
import argparse
from pathlib import Path
from typing import List, Tuple
import matplotlib.pyplot as plt
import prettypyplot as pplt  

def apply_pretty_style() -> None:
    for attr in ("set_style", "use_style", "apply_style"):
        fn = getattr(pplt, attr, None)
        if callable(fn):
            try:
                fn()
                return
            except TypeError:
                try:
                    fn("default")
                    return
                except Exception:
                    pass
            except Exception:
                pass

    # Enforce white backgrounds regardless of external style/theme.
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["savefig.edgecolor"] = "white"


def load_mc_energy(path: Path) -> Tuple[List[int], List[float]]:
    """Load accepted step and energy from dock_energy_mc.dat."""
    steps: List[int] = []
    energies: List[float] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                steps.append(int(parts[0]))
                energies.append(float(parts[1]))
            except ValueError:
                continue
    return steps, energies


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot MC docking energy trajectory."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="dock_energy_mc.dat",
        help="Input data file (default: dock_energy_mc.dat)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="dock_energy_mc.png",
        help="Output image file (default: dock_energy_mc.png)",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Optional max x-axis value (default: show full range)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    steps, energies = load_mc_energy(input_path)

    apply_pretty_style()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(steps, energies, color="#1f77b4", linewidth=1.5)
    ax.scatter(steps, energies, s=10, color="#1f77b4", alpha=0.6)
    ax.set_title("Monte Carlo Docking Energy Trajectory")
    ax.set_xlabel("Accepted MC Step")
    ax.set_ylabel("Complex Energy")
    ax.grid(True, alpha=0.25, linestyle="--", linewidth=0.8)
    if args.xmax is not None:
        ax.set_xlim(0, args.xmax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, facecolor="white", edgecolor="white")

    print(f"Saved figure: {output_path}")
    print(f"Loaded points: {len(steps)}")
    print(f"Min energy: {min(energies):.5f}")
    print(f"Final energy: {energies[-1]:.5f}")
if __name__ == "__main__":
    main()
