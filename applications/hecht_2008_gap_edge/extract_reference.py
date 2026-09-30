"""Decode selected published analytic/NRG spectral curves from Grace EPS vectors."""

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
SOURCE = "https://export.arxiv.org/e-print/0803.1251v3"


def figure_points(text, figure):
    paths = re.findall(r"(?m)^n\n((?:[-\d.]+ [-\d.]+ [ml]\n)+)s", text)
    if figure == 6:
        selections = ((3, "analytic", 1e-4, 0., "a"), (4, "NRG", 1e-4, 0., "a"),
                      (138, "analytic", 1e-4, 0., "b"), (139, "analytic", 1e-4, .016, "b"),
                      (140, "analytic", 1e-4, .032, "b"), (141, "NRG", 1e-4, .032, "b"),
                      (142, "analytic", 1e-4, .064, "b"))
    else:
        selections = ((2, "NRG", 3e-5, -.3, "c"), (5, "NRG", .03, -.3, "c"))
    rows = []
    for index, method, delta, epsilon, panel in selections:
        vertices = np.array([list(map(float, line.split()[:2])) for line in paths[index].splitlines()])
        if len(vertices) < 100:
            raise ValueError("unexpected published curve structure")
        for x, y in vertices[::5]:
            if figure == 7:
                logx = -8+(x-.3175)*8/(.9874-.3175)
                logy = -1+(y-.2180)*3/(.7306-.2180)
            elif panel == "a":
                logx = -15+(x-.2819)*10/(.9415-.2819)
                logy = (y-.3999)*2/(.6498-.3999)
            else:
                logx = -10+(x-1.8301)*5/(2.2153-1.8301)
                logy = -2+(y-.1681)*3/(.7281-.1681)
            rows.append(dict(figure=figure, panel=panel, method=method, delta_over_D=delta,
                             epsilon_over_D=epsilon, omega_prime_over_D=10**logx,
                             pi_gamma_A=10**logy, eps_x=x, eps_y=y))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    args = parser.parse_args()
    source = fetch(SOURCE)
    rows, members = [], {}
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        for number, name in ((6, "figure_6ab.eps"), (7, "figure_7cd.eps")):
            eps = archive.extractfile(name).read()
            rows.extend(figure_points(eps.decode(), number))
            members[name] = hashlib.sha256(eps).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"spectra.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1088/0953-8984/20/27/275213", arxiv_version="0803.1251v3",
        archive_url=SOURCE, archive_sha256=hashlib.sha256(source).hexdigest(), members_sha256=members,
        reference_type="Original log-log vector curves, every fifth vertex; not raw NRG spectra",
        coordinate_half_step=.00005, estimated_relative_coordinate_uncertainty=.003,
        notes=["Fig. 6(a): Gamma=0.008 D, epsilon=0; Fig. 6(b): epsilon/Gamma=0,2,4,8.",
               "Fig. 6(b) thick NRG curve is epsilon/Gamma=4 (dotted analytic curve).",
               "Fig. 7(c): U=0.6 D, Gamma=0.049 D, epsilon=-U/2; selected gaps 3e-5 and 0.03 D.",
               "Eq. (30) prints an extra Gamma in the squared denominator term. Use Eq. (28), the normal-state limit and the original curves to check the correction.",
               "NRG broadening and truncation errors are additional to vector-coordinate precision."]))
    print(f"Extracted {len(rows)} spectral vertices")


if __name__ == "__main__":
    main()
