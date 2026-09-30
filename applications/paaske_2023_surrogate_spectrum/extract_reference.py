"""Extract the published Fig. 3 vector paths and the authors' bath coefficients.

This optional maintenance command needs PyMuPDF: python -m pip install pymupdf.
Numerical runs and regression tests use the committed extracted data offline.
"""

import argparse
import csv
import hashlib
from io import BytesIO
from pathlib import Path
import tarfile
from urllib.parse import quote

import numpy as np

from applications._reference import fetch
from applications._common import write_json
from qdjj_solver import DiscreteBath

CASE = Path(__file__).resolve().parent
SOURCE = "https://export.arxiv.org/e-print/2307.11646v3"
COMMIT = "1ef00d1018a61be5a2b892ce3b9d2082288e81ec"
REPO = f"https://raw.githubusercontent.com/virgilb91/surrogate_models/{COMMIT}/"


def figure_points(drawings):
    """Calibrate original PDF path vertices from the vector tick positions."""
    rows = []
    for index in range(2, 16):
        path = drawings[index]
        if len(path["items"]) < 10 or any(item[0] != "l" for item in path["items"]):
            raise ValueError("published figure path structure changed")
        nrg = index >= 14
        label = "NRG" if nrg else f"L{(index-2)//2+1}"
        branch = "Dg" if index % 2 else "S"
        vertices = [path["items"][0][1]]+[item[2] for item in path["items"]]
        for point in vertices:
            if not 37.73 <= point.x <= 258.814 or not 20.366 <= point.y <= 123.169:
                continue
            gamma = 10**((point.x-89.5)/(219.594-89.5))
            excitation = (123.168-point.y)/((123.168-20.367)/2)
            rows.append(dict(method=label, branch=branch, gamma=float(gamma),
                             excitation=float(excitation), pdf_x=point.x, pdf_y=point.y))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    import pymupdf
    source = fetch(SOURCE)
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        pdf = archive.extractfile("fig_3.pdf").read()
    with pymupdf.open(stream=pdf, filetype="pdf") as document:
        page = document[0]
        if args.inspect:
            print("page", page.rect)
            print("text", page.get_text())
            print("words", page.get_text("words"))
            for i, drawing in enumerate(page.get_drawings()):
                if len(drawing["items"]) > 5:
                    print(i, "color", drawing["color"], "width", drawing["width"],
                          "dashes", drawing["dashes"], "rect", drawing["rect"],
                          "items", len(drawing["items"]), "first", drawing["items"][:2])
                elif drawing["color"] == (0., 0., 0.):
                    print("axis", i, drawing["items"])
            print("source sha256", hashlib.sha256(source).hexdigest())
            return
        rows = figure_points(page.get_drawings())
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure3.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    baths, files = {}, {}
    for levels in range(1, 9):
        arrays = []
        for prefix in ("xi", "gamma"):
            name = f"{prefix}_L={levels}_B=10_wc=10.dat"
            data = fetch(REPO+quote(name, safe=""))
            files[name] = hashlib.sha256(data).hexdigest()
            arrays.append(np.atleast_1d(np.loadtxt(BytesIO(data))))
        xi, residue = arrays
        order = np.argsort(xi)
        if xi.shape != (levels,) or residue.shape != xi.shape:
            raise ValueError("unexpected published bath layout")
        bath = DiscreteBath(xi[order], np.pi/20*residue[order], delta=1., bandwidth=10.,
                            metadata=dict(kind="published-surrogate", levels=levels, frequency_cutoff=10.,
                                          source_commit=COMMIT))
        baths[str(levels)] = bath.record()
    write_json(args.output/"baths.json", baths)
    (args.output/"LICENSE-bath-data.txt").write_text(fetch(REPO+"LICENSE").decode(), encoding="utf-8")
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PhysRevB.108.L220506", arxiv_version="2307.11646v3",
        source_archive=SOURCE, source_sha256=hashlib.sha256(source).hexdigest(),
        member="fig_3.pdf", member_sha256=hashlib.sha256(pdf).hexdigest(),
        calibration=dict(gamma1_x=89.5, gamma10_x=219.594, energy0_y=123.168, energy2_y=20.367),
        reference_type="Vector path coordinates, not the authors' raw NRG table",
        extraction_energy_uncertainty=0.0005, extraction_gamma_relative_uncertainty=0.0001,
        bath_repository=REPO, bath_files_sha256=files, bath_license="MIT",
        nrg=dict(Lambda=2, kept_states=500, z_averaging=False, half_bandwidth_over_delta=10),
        notes=["Dg is the coupled gerade doublet, which can lie above the continuum edge in a small surrogate.",
               "Colored lines are published surrogate-model calculations; black lines are independent NRG.",
               "Interpolation between sparse source vertices introduces additional comparison error."]))
    print(f"Extracted {len(rows)} figure vertices and eight published baths")


if __name__ == "__main__":
    main()
