"""Recover Fig. 1(a) NRG marker centers from the original MATLAB EPS."""

import argparse
import csv
import hashlib
from io import BytesIO
from pathlib import Path
import re
import tarfile

import numpy as np

from applications._common import write_json
from applications._reference import fetch

CASE = Path(__file__).resolve().parent
SOURCE = "https://export.arxiv.org/e-print/1610.08366v3"


def figure_points(text):
    paths = re.findall(r"\bN\s+((?:[-\d.]+ [-\d.]+ [ML]\s+)+)(cp\s+)?S", text)
    groups, points = [], []
    for path, closed in paths:
        vertices = np.array([list(map(float, m)) for m in re.findall(r"([-\d.]+) ([-\d.]+) [ML]", path)])
        if not np.all((vertices[:, 0] >= 169) & (vertices[:, 0] <= 447)):
            continue
        if closed and len(vertices) == 10 and np.isclose(np.ptp(vertices[:, 0]), 8.):
            # Ten-sided circular markers of radius four; export rounded to half-pixels.
            points.append((vertices[0, 0]-4., vertices[0, 1]))
        elif len(vertices) > 20 and points:
            groups.append(points)
            points = []
    if len(groups) != 9:
        raise ValueError(f"expected nine Fig. 1(a) curves, found {len(groups)}")
    rows = []
    for u, group in zip((2., 2.5, 2.8, 3., 3.2, 4., 5., 6., 7.), groups, strict=True):
        for x, y in sorted(group):
            rows.append(dict(u_meV=u, tilde_epsilon=(x-169)/278,
                             phi_over_pi=float(np.clip((338-y)/307, 0, 1)), eps_x=x, eps_y=y))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    args = parser.parse_args()
    source = fetch(SOURCE)
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        eps = archive.extractfile("fig_phase-boundary.eps").read()
    rows = figure_points(eps.decode())
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure1.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PhysRevB.95.195114", arxiv_version="1610.08366v3",
        archive_url=SOURCE, archive_sha256=hashlib.sha256(source).hexdigest(),
        member="fig_phase-boundary.eps", member_sha256=hashlib.sha256(eps).hexdigest(),
        reference_type="Original NRG marker centers, not interpolation curves or raw NRG output",
        calibration=dict(x_zero=169, x_one=447, y_zero=338, y_pi=31),
        tilde_epsilon_uncertainty=.0018, phi_over_pi_uncertainty=.0017,
        notes=["Coordinates before the global 0.75 scale and vertical reflection.",
               "Circle centers are exported on a half-pixel grid; endpoints clipped within rounding error.",
               "Curves ordered by increasing U as specified in the caption.",
               "Chi here is the squared phasor magnitude, unlike the Zalom 2024 application.",
               "No numerical bandwidth is specified for Fig. 1; finite-band sensitivity is measured separately."]))
    print(f"Extracted {len(rows)} NRG markers")


if __name__ == "__main__":
    main()
