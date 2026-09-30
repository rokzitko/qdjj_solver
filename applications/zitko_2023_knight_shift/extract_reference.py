"""Extract the authors' NRG numbers from their CC BY 4.0 Zenodo deposit.

Only the small resulting tables are retained. The two ZIP downloads total 88 MB.
No archived code or notebooks are executed. Run from the repository root:
python -m applications.zitko_2023_knight_shift.extract_reference
"""

import argparse
import csv
import hashlib
from io import BytesIO
from pathlib import Path
import urllib.request
import zipfile

import numpy as np

from applications._common import write_json

CASE = Path(__file__).resolve().parent
BASE = "https://zenodo.org/api/records/7951006/files/"
FILES = {"fig9.zip": "71b71aacf69ec161cf3044ec5a928f40",
         "fig10.zip": "728a54b21c9a3fecb5fd8a9d22a311f0"}


def grace_sets(text):
    """Grace numerical sets are terminated by &, independently of plot styling."""
    sets, rows = [], []
    for line in text.splitlines():
        line = line.strip()
        if line == "&":
            if rows:
                sets.append(np.array(rows))
                rows = []
        elif line and not line.startswith(("@", "#")):
            rows.append([float(x) for x in line.split()])
    if rows:
        sets.append(np.array(rows))
    return sets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    archives, provenance = {}, {}
    for name, checksum in FILES.items():
        with urllib.request.urlopen(BASE+name+"/content", timeout=180) as response:
            data = response.read()
        if hashlib.md5(data).hexdigest() != checksum:
            raise ValueError(f"archive checksum mismatch: {name}")
        archives[name] = zipfile.ZipFile(BytesIO(data))
        provenance[name] = dict(url=BASE+name+"/content", md5=checksum,
                                sha256=hashlib.sha256(data).hexdigest())
    member9 = "fig9/figphi1a.agr"
    member10 = "fig10/U=0.1-Gamma=1e-2-Lambda=4/corr.dat"
    sets = grace_sets(archives["fig9.zip"].read(member9).decode())
    text10 = archives["fig10.zip"].read(member10).decode()
    if args.inspect:
        for i, data in enumerate(sets):
            print("set", i, "shape", data.shape, "first", data[:3], "last", data[-2:])
        print("Fig10 corr.dat:", text10[:3000])
        return
    if len(sets) != 2 or not np.allclose(sets[0][:, 0], sets[1][:, 0]):
        raise ValueError("unexpected Fig. 9 data layout")
    phase = np.loadtxt(BytesIO(text10.encode()))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"figure9.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["gamma_over_u", "kappa_phi0", "kappa_phipi"])
        writer.writerows(zip(sets[0][:, 0], sets[0][:, 1], sets[1][:, 1], strict=True))
    with (args.output/"figure10.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["nominal_phi_over_pi", "phi_over_pi", "kappa"])
        writer.writerows((p, p+0.0001, k) for p, k in phase[:, :2])
    write_json(args.output/"provenance.json", dict(
        doi="10.5281/zenodo.7951006", license="CC-BY-4.0",
        creators=["Luka Pavešić", "Rok Žitko"], archives=provenance,
        members={"figure9.csv": member9, "figure10.csv": member10},
        notes=["Gamma is total physical hybridization.",
               "Fig. 9 uses Lambda=8 and averages z=1,0.5; E_Z/Delta=0.01.",
               "Fig. 10 uses Lambda=4,z=1; E_Z/Delta=0.001.",
               "corr.dat contains nominal phases; actual phi/pi is nominal+0.0001 as in the archived sweep."]))


if __name__ == "__main__":
    main()
