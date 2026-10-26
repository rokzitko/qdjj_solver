"""Shared models, data records and figure style for the physics manual.

Each ``fig_*.py`` script has two subcommands. ``compute`` runs the solver and
writes a JSON record to ``docs/manual/data`` (tracked, with provenance);
``plot`` turns that record into a PDF figure in ``docs/manual/figures``
(an ignored build artifact). Plotting never runs the solver.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA, FIGURES = ROOT/"data", ROOT/"figures"

# Reference junction of the manual: energies in units of Delta=1.
REFERENCE = json.loads((ROOT.parents[1]/"examples"/"01_subgap_states"/"input"/"reference.json").read_text())
BATH = {k: v for k, v in REFERENCE["bath"].items() if k != "kind"}
JUNCTION = {k: v for k, v in REFERENCE["model"].items() if k != "kind"}


def bath(**changes):
    from qdjj_solver import cosh_grid
    return cosh_grid(**{**BATH, **changes})


def junction(**changes):
    """Keyword arguments of ``reference_model`` for the manual's junction."""
    base = dict(JUNCTION)
    pairs = changes.pop("pairs", BATH["pairs"])
    bandwidth = changes.pop("bandwidth", BATH["bandwidth"])
    base.update(changes)
    base["bath"] = bath(pairs=pairs, bandwidth=bandwidth)
    return base


def with_omega(periodic, omega):
    """The same harmonics at another frequency (exact for H(p0 + A cos(omega t)))."""
    from qdjj_solver import PeriodicHamiltonian
    return PeriodicHamiltonian(omega, periodic.harmonics, periodic.nimp, periodic.spins, periodic.eta_labels,
                               periodic.reference, periodic.observables, periodic.metadata)


def basis_size(energies, minimum):
    """Smallest static-basis size K >= minimum that does not split a degenerate multiplet."""
    gaps = np.flatnonzero(np.diff(energies) > 1e-7) + 1
    gaps = gaps[gaps >= minimum]
    return int(gaps[0]) if len(gaps) else len(energies)


def provenance(script):
    from qdjj_solver import __version__
    from qdjj_solver.common.io import numerical_environment
    try:
        commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                                    capture_output=True, text=True, check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    return dict(script=Path(script).name, qdjj_solver=__version__, git_commit=commit, src_modified=dirty,
                environment=numerical_environment(), created=time.strftime("%Y-%m-%dT%H:%M:%S%z"))


def save(name, payload, script, compact=False):
    """Write a record with provenance; ``compact`` drops indentation for large point clouds."""
    from qdjj_solver.common.io import json_value, write_json
    DATA.mkdir(parents=True, exist_ok=True)
    record = dict(provenance=provenance(script), **payload)
    if compact:
        text = json.dumps(json_value(record), sort_keys=True, separators=(",", ":"), allow_nan=False)
        (DATA/f"{name}.json").write_text(text + "\n", encoding="utf-8")
    else:
        write_json(DATA/f"{name}.json", record)
    print(f"wrote {DATA/f'{name}.json'}")


def load(name):
    return json.loads((DATA/f"{name}.json").read_text(encoding="utf-8"))


def main(compute, plot):
    parser = argparse.ArgumentParser(description=sys.modules["__main__"].__doc__)
    parser.add_argument("action", choices=("compute", "plot", "all"))
    action = parser.parse_args().action
    if action in ("compute", "all"):
        start = time.perf_counter()
        compute()
        print(f"compute: {time.perf_counter()-start:.1f} s")
    if action in ("plot", "all"):
        plot()


# Figure style: light print surface, validated categorical slots 1-3, one-hue sequential ramp.
# Notation follows the text: \varepsilon is a quasienergy and the lunate ϵ (U+03F5) the detuning, as \eps
# in the manual. Mathtext draws \epsilon like \varepsilon, so labels write the detuning as ϵ directly.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, SECONDARY, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
WIDTH = 5.9  # inches; text width of the manual


def style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.edgecolor": AXIS, "axes.linewidth": .6, "axes.labelcolor": INK, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": SECONDARY, "ytick.labelcolor": SECONDARY,
        "xtick.major.width": .6, "ytick.major.width": .6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": .4, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
        "lines.solid_capstyle": "round", "legend.frameon": False, "legend.labelcolor": INK,
        "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": .02, "pdf.fonttype": 42,
    })


def figure(rows=1, columns=1, height=2.3, width=WIDTH, **options):
    import matplotlib.pyplot as plt
    style()
    return plt.subplots(rows, columns, figsize=(width, height), constrained_layout=True, **options)


def panel_label(axis, label):
    axis.text(-.02, 1.02, f"({label})", transform=axis.transAxes, ha="right", va="bottom",
              fontsize=9, fontweight="bold", color=INK)


def finish(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES/f"{name}.pdf")
    print(f"wrote {FIGURES/f'{name}.pdf'}")


def ramp_color(value):
    """Sequential blue for a weight in [0, 1]; small weights recede toward the surface."""
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("manual-blue", RAMP)
    return cmap(np.clip(value, 0., 1.))


def weight_colormap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("manual-blue", ["#f0efec"] + RAMP)


def weighted_line(axis, x, y, weight, color, width=1.8, floor=.04, zorder=2):
    """A line whose opacity and width follow a weight in [0, 1] point by point."""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import to_rgba
    points = np.column_stack([x, y])
    segments = np.stack([points[:-1], points[1:]], axis=1)
    middle = np.clip((np.asarray(weight)[:-1] + np.asarray(weight)[1:])/2, 0., 1.)
    keep = middle > floor
    rgba = np.array([to_rgba(color, alpha=.15 + .85*np.sqrt(w)) for w in middle[keep]])
    collection = LineCollection(segments[keep], colors=rgba, linewidths=.4 + width*np.sqrt(middle[keep]),
                                capstyle="round", zorder=zorder)
    axis.add_collection(collection)
    return collection
