"""Save completed phase-scan points and plot from their scalar records.

python examples/handbook/scan.py --output runs/handbook/scan
python examples/handbook/scan.py --output runs/handbook/scan --plot-only
"""

import argparse
import json
from pathlib import Path

import numpy as np

from qdjj_solver import Sector, save_result, solve
from qdjj_solver.common.io import model_from_config, write_json

INPUT = Path(__file__).resolve().parents[1]/"01_subgap_states"/"input"/"reference.json"


def compute(config, output, points):
    records = []
    # BEGIN handbook: scan
    for index, phi in enumerate(np.linspace(0., np.pi, points)):
        point = {**config, "model": {**config["model"], "phi": float(phi)}}
        h = model_from_config(point, default_compress=False)
        row = dict(index=index, phi=float(phi), status="ok")
        try:
            even = solve(h, cutoff=point["cutoff"], sector=Sector(0, 0))
            odd = solve(h, cutoff=point["cutoff"], sector=Sector(1, 1))
            save_result(even, output/f"{index:03d}-even.json")
            save_result(odd, output/f"{index:03d}-odd.json")
            row.update(singlet=float(even.energies[0]), doublet=float(odd.energies[0]),
                       current_singlet=float(np.real(even.observables["phase_derivative"][0])),
                       current_doublet=float(np.real(odd.observables["phase_derivative"][0])),
                       residual=max(float(even.residuals[0]), float(odd.residuals[0])))
        except (ValueError, RuntimeError) as error:
            row.update(status="failed", error=str(error))
        records.append(row)
        write_json(output/"scan.json", dict(input=config, points=records))
    # END handbook: scan
    if any(row["status"] != "ok" for row in records):
        raise RuntimeError("scan has failed points; inspect scan.json before plotting")


def plot(output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = json.loads((output/"scan.json").read_text())["points"]
    if any(row["status"] != "ok" for row in rows):
        raise ValueError("cannot plot a complete scan with failed points")
    x = np.array([r["phi"] for r in rows])/np.pi
    gap = np.array([r["doublet"]-r["singlet"] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.5), layout="constrained")
    axes[0].plot(x, gap, "o-", color="#2a78d6")
    axes[0].axhline(0., color="gray", lw=.6)
    axes[0].set_ylabel(r"$(E_D-E_S)/\Delta$")
    for key, label, color in (("current_singlet", "S", "#2a78d6"), ("current_doublet", "D", "#eb6834")):
        axes[1].plot(x, [r[key] for r in rows], "o-", label=label, color=color)
    axes[1].set_ylabel(r"$I/(2e\Delta/\hbar)$")
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel(r"$\phi/\pi$")
    fig.savefig(output/"scan.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=Path("results/handbook-scan"))
    parser.add_argument("--points", type=int, default=9)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if args.points < 2:
        parser.error("--points must be at least two")
    if not args.plot_only:
        compute(json.loads(args.input.read_text()), args.output, args.points)
    plot(args.output)


if __name__ == "__main__":
    main()
