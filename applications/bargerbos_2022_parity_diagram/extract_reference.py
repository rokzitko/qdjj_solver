"""Extract the deposited processed NRG gate--phase map (4TU, CC BY 4.0).

Only selected ZIP members are transferred using verified HTTP Range requests.
A restricted numeric-array decoder reads the legacy pickle. No notebooks or
deposited programs are executed.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

from applications._common import write_json
from applications._reference import RemoteZip, numeric_pickle

CASE = Path(__file__).resolve().parent
URL = "https://data.4tu.nl/file/bd5acfda-aaee-49e8-8d92-9dc670c7374d/40b032ed-f590-4cac-849d-1b7c2eaae095"
ETAG = "487aa316f331ebc600694e01259c249a"
MEMBER = "Processed_data/ABS_map_xi_vs_phi_Gamma_0.20_xi_01_U05_Asym1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CASE/"reference")
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    remote = RemoteZip(URL, ETAG)
    with zipfile.ZipFile(remote) as archive:
        if args.inspect:
            print("Members:", [(i.filename, i.file_size) for i in archive.infolist()
                               if "ABS_map" in i.filename or i.filename.endswith("ipynb")])
        data = archive.read(MEMBER)
        gate, phase, gamma, energy_1, energy_2 = numeric_pickle(data)
        if args.inspect:
            for name, a in zip(("gate", "phase", "gamma", "energy_1", "energy_2"),
                               (gate, phase, gamma, energy_1, energy_2), strict=True):
                a = np.asarray(a)
                print(name, a.shape, "range", a.min(), a.max(), "first", a.reshape(-1)[:6])
            notebook = json.loads(archive.read("4-plot_figures_device_A.ipynb"))
            for cell in notebook["cells"]:
                source = "".join(cell.get("source", []))
                if "ABS_map_xi_vs_phi" in source or "'Fig_S1.pdf'" in source:
                    print("Notebook source (read only):", source)
            print("Bytes transferred:", remote.transferred)
            return
    if not np.isclose(gamma, .2) or energy_1.shape != (len(gate), len(phase)) or energy_2.shape != energy_1.shape:
        raise ValueError("unexpected reference parameters or array layout")
    rows = []
    for j, p in enumerate(phase):
        if not 0 <= p <= 1:
            continue
        gap = energy_1[:, j]-energy_2[:, j]
        crossings = [i for i in range(len(gate)-1) if gate[i] >= 0 and gap[i]*gap[i+1] < 0]
        if len(crossings) != 1:
            raise ValueError("expected exactly one positive-gate crossing")
        i = crossings[0]
        root = gate[i]-gap[i]*(gate[i+1]-gate[i])/(gap[i+1]-gap[i])
        rows.append(dict(phi_over_pi=float(p), detuning_critical_over_u=float(root),
                         gate_lower=float(gate[i]), gate_upper=float(gate[i+1]),
                         gap_lower=float(gap[i]), gap_upper=float(gap[i+1])))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/"boundary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output/"provenance.json", dict(
        paper_doi="10.1103/PRXQuantum.3.030311", dataset_doi="10.4121/19102769.v1",
        creators=["Arno Bargerbos", "Marta Pita Vidal"], license="CC-BY-4.0",
        archive_url=URL, archive_published_md5_and_verified_etag=ETAG,
        member=MEMBER, member_sha256=hashlib.sha256(data).hexdigest(),
        retrieval_bytes=remote.transferred,
        source="Processed/interpolated NRG grid used for supplementary Fig. S1(c)",
        operation="Linear zero crossing between adjacent positive xi/U grid points",
        gate_spacing_over_u=.01, phase_spacing_over_pi=.04,
        model=dict(u_over_delta=5, gamma_total_over_delta=1, half_bandwidth_over_delta=10),
        nrg=dict(Lambda=8, maximum_spin_multiplets=3000, z_averaging=False),
        notes=["Each member's ZIP CRC is verified; the entire archive is not downloaded or independently hashed.",
               "Zero crossings use the stored A1-A2 difference without assuming its sector order.",
               "Grid interpolation and NRG discretization errors exceed floating-point precision."]))
    print(f"Extracted {len(rows)} phase-boundary points; transferred {remote.transferred} bytes")


if __name__ == "__main__":
    main()
