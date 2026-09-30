"""Extract original vector coordinates of Figs. 8 and 17 (optional PyMuPDF tool)."""

import argparse
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
import tarfile

from applications._reference import fetch

CASE = Path(__file__).resolve().parent
SOURCE = "https://export.arxiv.org/e-print/2508.18465v2"
FIGURE17_CALIBRATION = ((196.37380981445312, 43.50616455078125, 1.),
                        (308.9130554199219, 234.24652099609375, .1),
                        (427.9520568847656, 503.26849365234375, -.5))


def color_key(color):
    return tuple(round(x, 3) for x in color) if color is not None else None


def line_vertices(drawing):
    points = []
    for item in drawing["items"]:
        if item[0] == "l":
            for point in item[1:]:
                xy = (point.x, point.y)
                if not points or xy != points[-1]:
                    points.append(xy)
    return points


def figure8(page):
    drawings = page.get_drawings()
    frames = [(i, d["rect"]) for i, d in enumerate(drawings) if d["type"] == "f"
              and color_key(d["fill"]) == (1., 1., 1.)
              and 299 < d["rect"].width < 300 and 250 < d["rect"].height < 251]
    if len(frames) != 6:
        raise ValueError("expected six Fig. 8 frames")
    methods = ("GAL", "ChE-F1", "ChE-W2", "ChE-F2", "ChE-W4", "ChE-F4")
    colors = {(0., 0., 0.): 2., (1., 0., 0.): 4., (0., 0., 1.): 8.}
    rows, markers = [], {}
    for panel, (start, frame) in enumerate(frames):
        stop = frames[panel+1][0] if panel+1 < len(frames) else len(drawings)
        # Axes: -0.2 at the bottom, +0.7 at the top; phi/pi from zero to one.
        def physical(x, y, frame=frame):
            return (x-frame.x0)/frame.width, -.2+.9*(frame.y1-y)/frame.height
        for drawing in drawings[start+1:stop]:
            color = color_key(drawing["color"])
            if color not in colors:
                continue
            if drawing["type"] == "s" and len(drawing["items"]) >= 10:
                for x, y in line_vertices(drawing):
                    phase, current = physical(x, y)
                    if -1e-6 <= phase <= 1+1e-6:
                        rows.append(dict(method=methods[panel], u=colors[color],
                                         phi_over_pi=min(1., max(0., phase)), current=current,
                                         pdf_x=x, pdf_y=y))
            # Panel (a) contains two identical marker layers. Exclude legend
            # markers, drawn AFTER the axis ticks and the GAL polylines.
            if panel == 0 and drawing["type"] == "fs" and len(drawing["items"]) == 1:
                break
            if panel == 0 and drawing["type"] == "fs" and len(drawing["items"]) in (6, 8):
                rect = drawing["rect"]
                x, y = (rect.x0+rect.x1)/2, (rect.y0+rect.y1)/2
                phase, current = physical(x, y)
                key = (colors[color], round(phase, 6))
                markers[key] = dict(method="NRG", u=colors[color], phi_over_pi=key[1],
                                    current=current, pdf_x=x, pdf_y=y)
    # Collect GAL lines separately: the marker-layer termination above is
    # intentional and prevents the inset legend from becoming fake NRG data.
    start, frame = frames[0]
    for drawing in drawings[start+1:frames[1][0]]:
        color = color_key(drawing["color"])
        if drawing["type"] == "s" and color in colors and len(drawing["items"]) >= 10:
            for x, y in line_vertices(drawing):
                rows.append(dict(method="GAL", u=colors[color], phi_over_pi=(x-frame.x0)/frame.width,
                                 current=-.2+.9*(frame.y1-y)/frame.height, pdf_x=x, pdf_y=y))
    rows.extend(markers.values())
    if len(markers) != 72:
        raise ValueError(f"unexpected number of distinct NRG markers: {len(markers)}")
    return sorted(rows, key=lambda r: (r["method"], r["u"], r["phi_over_pi"]))


def figure17(page):
    drawings = page.get_drawings()
    frames = [(i, d["rect"]) for i, d in enumerate(drawings) if d["type"] == "f"
              and color_key(d["fill"]) == (1., 1., 1.) and 557 < d["rect"].width < 559]
    if len(frames) != 3:
        raise ValueError("expected three Fig. 17 frames")
    colors = {(1., .647, 0.): "ChE-W2", (.216, .494, .722): "ChE-W4",
              (0., .502, 0.): "ChE-W6", (1., 0., 0.): "ChE-W8"}
    # Independent axis-tick calibration, PDF coordinates (origin at top left).
    rows = []
    for panel, (start, frame) in enumerate(frames):
        stop = frames[panel+1][0] if panel+1 < len(frames) else len(drawings)
        zero, tick, value = FIGURE17_CALIBRATION[panel]
        branches = {}
        for drawing in drawings[start+1:stop]:
            color = color_key(drawing["color"])
            if drawing["type"] == "fs" and color == (0., 0., 0.) and len(drawing["items"]) == 8:
                rect = drawing["rect"]
                points = [((rect.x0+rect.x1)/2, (rect.y0+rect.y1)/2)]
                method, branch = "NRG", -1
            elif drawing["type"] == "s" and color in colors and len(drawing["items"]) > 10:
                points = line_vertices(drawing)
                method = colors[color]
                branch = branches.get(method, 0)
                branches[method] = branch+1
            else:
                continue
            for x, y in points:
                if not frame.x0-.01 <= x <= frame.x1+.01 or not frame.y0-.01 <= y <= frame.y1+.01:
                    continue
                rows.append(dict(method=method, observable=("excitation", "pairing", "spin_correlation")[panel],
                                 branch=branch, u=round((x-frame.x0)/frame.width*10, 6),
                                 value=(y-zero)*value/(tick-zero), pdf_x=x, pdf_y=y))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    parser.add_argument("--preview", type=Path, help="optional directory for original figure previews")
    args = parser.parse_args()
    import pymupdf
    source = fetch(SOURCE)
    members, tables, calibrations = {}, {}, {}
    with tarfile.open(fileobj=BytesIO(source)) as archive:
        for name, extractor in (("fig_08.pdf", figure8), ("fig_17.pdf", figure17)):
            data = archive.extractfile(name).read()
            members[name] = hashlib.sha256(data).hexdigest()
            document = pymupdf.open(stream=data, filetype="pdf")
            tables[name] = extractor(document[0])
            width = (299, 300) if name == "fig_08.pdf" else (557, 559)
            frames = [list(d["rect"]) for d in document[0].get_drawings() if d["type"] == "f"
                      and color_key(d["fill"]) == (1., 1., 1.) and width[0] < d["rect"].width < width[1]]
            calibrations[name] = dict(pdf_frames=frames, x_bounds=[0., 1. if name == "fig_08.pdf" else 10.])
            if name == "fig_08.pdf":
                calibrations[name].update(y_bounds=[-.2, .7], panel_order=["GAL", "ChE-F1", "ChE-W2", "ChE-F2", "ChE-W4", "ChE-F4"])
            else:
                calibrations[name].update(y_zero_tick_value=FIGURE17_CALIBRATION,
                                           panel_order=["excitation", "pairing", "spin_correlation"])
            if args.preview:
                args.preview.mkdir(parents=True, exist_ok=True)
                document[0].get_pixmap(matrix=pymupdf.Matrix(1.5, 1.5)).save(args.preview/(name+".png"))
    args.output.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        with (args.output/("figure8.csv" if name == "fig_08.pdf" else "figure17.csv")).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    provenance = dict(paper_doi="10.1103/mxsl-fc96", arxiv_version="2508.18465v2",
                      archive_url=SOURCE, archive_sha256=hashlib.sha256(source).hexdigest(),
                      members_sha256=members, calibration=calibrations, license="CC BY 4.0 (arXiv version)",
                      reference_type="Original PDF marker centers and polyline vertices, not raw NRG output",
                      coordinate_precision_points=.001, figure8_phase_uncertainty=.000004,
                      figure8_current_uncertainty=.000004, figure17_u_uncertainty=.00002,
                      figure17_value_uncertainty=.00001,
                      notes=["Fig. 8 figure legend assigns black U=2, red U=4, blue U=8; the prose reverses the last two colors.",
                             "Fig. 8 NRG markers are taken from panel (a), with duplicate marker layers and legends removed.",
                             "Fig. 17 NRG multiplets have no spin labels in the figure; comparison uses one-to-one matching of displayed energy levels.",
                             "Curve simplification/interpolation error exceeds marker-coordinate rounding and is reported separately.",
                             "ChE-W denotes wide-band coefficients; ChE-F uses D/Delta=100.",
                             "NRG used D/Delta=100, Lambda=4 for two channels and Lambda=2 for one; kept-state and z-averaging settings not specified."])
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8")
    print({name: len(rows) for name, rows in tables.items()})


if __name__ == "__main__":
    main()
