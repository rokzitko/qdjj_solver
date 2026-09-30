"""Extract NRG marker centers from the arXiv vector version of Fig. 9(a).

This is a graphical extraction, not the authors' raw NRG data. The SVG is
calibrated by the frame corners: x=588..5039 -> phi/pi=0..1,
y=5188..6854.5 -> J/J0=-0.1..0.4, before its common display transform.
The two marker keys at x=902 are legend entries and are excluded explicitly.
"""

import argparse
import csv
import hashlib
from pathlib import Path
import re
import urllib.request
import xml.etree.ElementTree as ET

from applications._common import write_json

CASE = Path(__file__).resolve().parent
URL = "https://arxiv.org/html/2211.10312v2/fig_09.svg"


def extract(data):
    pattern = re.compile(r"M([\d.]+) ([\d.]+)V([\d.]+)$")
    points = set()
    for path in ET.fromstring(data).iter("{http://www.w3.org/2000/svg}path"):
        match = pattern.fullmatch(path.get("d", ""))
        if path.get("stroke-linecap") != "round" or match is None or match[2] != match[3]:
            continue
        x, y = float(match[1]), float(match[2])
        color = path.get("stroke")
        if color not in ("#0000ff", "#ff0000") or x == 902.:
            continue
        if 588 <= x <= 5039 and 5188 <= y <= 6854.5:
            points.add(("weak" if color == "#0000ff" else "strong", x, y))
    rows = [dict(curve=label, phi_over_pi=round((x-588)/4451, 3),
                 current=-.1+(y-5188)/1666.5*.5, svg_x=x, svg_y=y)
            for label, x, y in sorted(points)]
    if len(rows) != 69:
        raise ValueError(f"expected 69 distinct data markers, found {len(rows)}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    args = parser.parse_args()
    with urllib.request.urlopen(URL, timeout=60) as response:
        data = response.read()
    rows = extract(data)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure9a.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PhysRevB.107.115407", source=URL,
        source_sha256=hashlib.sha256(data).hexdigest(), figure="9(a)",
        kind="NRG marker centers extracted from vector figure",
        calibration=dict(x=[588, 5039], phase_over_pi=[0, 1],
                         y=[5188, 6854.5], current=[-.1, .4]),
        coordinate_rounding_half_step=.25,
        phase_over_pi_extraction_uncertainty=0.0001,
        current_extraction_uncertainty=0.00015,
        notes=["The plot quantizes coordinates; these are not raw NRG outputs.",
               "Identical duplicate marker removed; the two legend markers at x=902 excluded.",
               "Nominal phase grid recovered to 0.001 pi; SVG coordinates are retained.",
               "NRG calculation uncertainty is additional to graphical extraction uncertainty."]))
    print(f"Extracted {len(rows)} NRG points")


if __name__ == "__main__":
    main()
