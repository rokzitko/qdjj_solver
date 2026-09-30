"""Recover Fig. 3(a--c) current curves from the original Grace EPS source.

Coordinates are read as data, without executing PostScript. Published physical
parameters and units are kept exactly; no fitting of Gamma or current scale.
"""

import argparse
import csv
import hashlib
from io import BytesIO
from pathlib import Path
import re
import tarfile

from applications._common import write_json
from applications._reference import fetch

CASE = Path(__file__).resolve().parent
SOURCE = "https://export.arxiv.org/e-print/cond-mat/0312271v2"


def figure_points(text):
    rows, points, color = [], [], None
    number = r"(-?\d+\.\d+)"
    for line in text.splitlines():
        line = line.strip()
        match = re.fullmatch(r"Color(\d+) SC", line)
        if match:
            color = int(match[1])
        marker = re.fullmatch(rf"n {number} {number} 0\.0050 0\.0050 0 360 EARC s", line)
        if marker:
            x, y = map(float, marker.groups())
            if .15 <= x <= .4 and .9 <= y <= 1.1:
                rows.append(dict(delta_over_tk=10., phi_over_pi=(x-.275)/.125,
                                 current=(y-1)/10, eps_x=x, eps_y=y))
            elif .425 <= x <= .675 and .9 <= y <= 1.1:
                rows.append(dict(delta_over_tk=.1, phi_over_pi=(x-.55)/.125,
                                 current=(y-1)*10, eps_x=x, eps_y=y))
        if line == "n":
            points = []
        vertex = re.fullmatch(rf"{number} {number} ([ml])", line)
        if vertex:
            x, y = map(float, vertex.groups()[:2])
            if vertex[3] == "m":
                points = []
            points.append((x, y))
        if line == "s":
            if (len(points) >= 20 and color in (1, 2, 3, 4) and
                    all(.1499 <= x <= .6751 and .5499 <= y <= .8501 for x, y in points)):
                for x, y in points:
                    rows.append(dict(delta_over_tk={1: 1.6, 2: 1.8, 3: 2., 4: 2.2}[color],
                                     phi_over_pi=(x-.4125)/.2625, current=y-.7, eps_x=x, eps_y=y))
            points = []
    if {r["delta_over_tk"] for r in rows} != {.1, 1.6, 1.8, 2., 2.2, 10.}:
        raise ValueError("failed to identify all six published CPRs")
    return sorted(rows, key=lambda r: (r["delta_over_tk"], r["phi_over_pi"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    args = parser.parse_args()
    source = fetch(SOURCE)
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        eps = archive.extractfile("fig3.eps").read()
    rows = figure_points(eps.decode())
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure3.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PhysRevB.70.020502", arxiv_version="cond-mat/0312271v2",
        archive_url=SOURCE, archive_sha256=hashlib.sha256(source).hexdigest(),
        member="fig3.eps", member_sha256=hashlib.sha256(eps).hexdigest(),
        reference_type="Original vector-figure coordinates, not raw NRG output",
        calibration=dict(panel_a=dict(x_zero=.275, x_pi=.4, y_zero=1., current_per_y=.1),
                         panel_b=dict(x_zero=.55, x_pi=.675, y_zero=1., current_per_y=10.),
                         panel_c=dict(x_zero=.4125, x_pi=.675, y_zero=.7, current_per_y=1.)),
        coordinate_half_step=.00005, phase_over_pi_uncertainty=.0004,
        maximum_current_extraction_uncertainty=.0005,
        notes=["Panels (a,b) use circle centers; panel (c) uses colored polyline vertices.",
               "Current is quoted in e Delta/hbar. Physical Gamma is total coupling.",
               "Kondo scale follows paper Eq. (4), not an independently measured Kondo temperature.",
               "Original NRG results have additional discretization/truncation errors."]))
    print(f"Extracted {len(rows)} points from six current curves")


if __name__ == "__main__":
    main()
