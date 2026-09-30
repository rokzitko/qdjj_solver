"""Decode original Gnuplot Fig. 8 vectors without executing PostScript."""

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
SOURCE = "https://export.arxiv.org/e-print/1509.06959v2"


def figure_points(text):
    blocks = re.findall(r"% Begin plot #(\d+)(.*?)% End plot #\1", text, re.S)[:17]
    if len(blocks) != 17:
        raise ValueError("unexpected Fig. 8 main-panel structure")
    rows = []
    for index, ratio, method in ((1, 1., "FDC_L"), (3, 1., "DC_L"), (4, 1., "DC_minus_R"),
                                 (5, .5, "FDC_L"), (7, .5, "DC_L"), (8, .5, "DC_minus_R"),
                                 (9, .25, "FDC_L"), (11, .25, "DC_L"), (12, .25, "DC_minus_R"),
                                 (13, .5, "NRG"), (14, .5, "NRG")):
        body = blocks[index][1]
        points, x, y = [], 0., 0.
        for match in re.finditer(r"(-?\d+) (-?\d+) (M|V|CircleF)\b", body):
            a, b, op = match.groups()
            if op in ("M", "CircleF"):
                x, y = float(a), float(b)
            else:
                x, y = x+float(a), y+float(b)
            points.append((x, y))
        # The first FDC curve contains a legend before the actual data moveto.
        if index == 1:
            start = next(i for i, (x, _) in enumerate(points) if x == 1170)
            points = points[start:]
        for x, y in points:
            rows.append(dict(gap_ratio=ratio, method=method, detuning=-4+8*(x-1170)/5489,
                             current=(y-1110)/13995, eps_x=x, eps_y=y))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    args = parser.parse_args()
    source = fetch(SOURCE)
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        eps = archive.extractfile("J_LR.eps").read()
    rows = figure_points(eps.decode())
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure8.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PhysRevB.93.024523", arxiv_version="1509.06959v2",
        archive_url=SOURCE, archive_sha256=hashlib.sha256(source).hexdigest(),
        member="J_LR.eps", member_sha256=hashlib.sha256(eps).hexdigest(),
        reference_type="Original vector figure, not raw solver output",
        calibration=dict(x_minus4=1170, x_plus4=6659, y_zero=1110, y_point2=3909),
        coordinate_half_step=.5, gate_uncertainty=.00073, current_uncertainty=.000036,
        notes=["NRG circles trace the Delta_R/Delta_L=0.5 curve on both sides.",
               "FDC is plotted for negative detuning, DC for positive detuning.",
               "The paper sets e=hbar=1: J/Delta_L is in e Delta_L/hbar.",
               "Paper uses Phi_L-Phi_R; exchanging the phase orientation preserves the plotted positive CPR branch.",
               "The continuum formulas use a wide band; finite D is a separate convergence check."]))
    print(f"Extracted {len(rows)} Fig. 8 coordinates")


if __name__ == "__main__":
    main()
