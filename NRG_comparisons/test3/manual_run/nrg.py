"""Standalone SIAM/MAGFIELD/SPU1 decks and read-only final-shell extraction.

NRG Ljubljana: hamiltonian.m adds +B*spinz[d[]] to delta*n+U/2*(n-1)^2;
operators.m defines SZd and SZ2d as scalar operators; invar.hpp has
SZ(SSZ)=(SSZ-1)/2. No SU(2) reduced-matrix conversion applies here.
Only extract_result needs NumPy/h5py. physics.py is the shared local contract.
"""

from configparser import ConfigParser, Error as ConfigError
import math
from pathlib import Path
import re


FIXED = {
    "model": "SIAM", "variant": "MAGFIELD", "symtype": "SPU1", "band": "flat",
    "Ninit": "0", "lastall": "true", "discretization": "Z", "tri": "old",
    "options": "SCHUR", "diagratio": "1", "diag_mode": "serial", "diagth": "1",
    "strategy": "kept", "ops": "n_d n_d_ud SZd SZ2d", "finite": "true",
    "dm": "false", "fdm": "false", "broaden": "false", "savebins": "false",
    "calc0": "true", "h5raw": "true", "h5all": "false", "h5last": "true",
    "h5ops": "true", "h5vectors": "false", "h5U": "false",
    "removefiles": "true", "silent": "true",
}
NUMERICAL = {"lambda": "Lambda", "z": "z", "nmax": "Nmax", "keep": "keep",
             "keepmin": "keepmin", "keepenergy": "keepenergy"}
ENERGIES = {"bandwidth": "bandrescale", "u": "U", "gamma": "Gamma", "detuning": "delta"}


def render_param(physical, numerical, units="gap") -> str:
    """Render explicit settings; optional untruncated=True is for tiny controls.

    Energy inputs share a physical unit; the deck unit is gap or bandwidth.
    keepenergy uses Wilson-shell units. No Tmin or polarized bath is accepted.
    """
    if set(physical) != {*ENERGIES, "gap", "field", "geometry"}:
        raise ValueError("physical requires exactly geometry, gap, bandwidth, u, gamma, detuning, field")
    if not set(NUMERICAL) <= numerical.keys() or numerical.keys() - {*NUMERICAL, "untruncated"}:
        raise ValueError("unsupported or missing numerical settings")
    if units not in ("gap", "bandwidth") or physical["geometry"] != "single":
        raise ValueError("units must be gap or bandwidth; geometry must be single")
    for label, values, keys in (("physical", physical, {*ENERGIES, "gap", "field"}),
                                ("numerical", numerical, set(NUMERICAL))):
        for key in keys:
            value = values[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{label}.{key} must be finite and real")
    if min(physical["gap"], physical["bandwidth"]) <= 0 or min(physical["u"], physical["gamma"]) < 0:
        raise ValueError("gap/bandwidth must be positive; u/gamma must be nonnegative")
    if numerical["lambda"] < 1.8 or not 0 < numerical["z"] <= 1:
        raise ValueError("Lambda must be >= 1.8 and z in (0, 1]")
    for key, minimum in (("nmax", 1), ("keep", 2), ("keepmin", 0)):
        if type(numerical[key]) is not int or numerical[key] < minimum:
            raise ValueError(f"{key} must be an integer >= {minimum}")
    if numerical["nmax"] > 998 or numerical["keepmin"] > numerical["keep"]:
        raise ValueError("Nmax must be <= 998 and keepmin <= keep")
    if type(numerical.get("untruncated", False)) is not bool:
        raise ValueError("untruncated must be boolean")
    unit = physical[units]
    params = dict(FIXED)
    params.update({deck: format(physical[key] / unit, ".17g") for key, deck in ENERGIES.items()})
    params["bcsgap"] = format(physical["gap"] / physical["bandwidth"], ".17g")
    params.update({deck: format(numerical[key], ".17g") for key, deck in NUMERICAL.items()})
    # truncation.hpp uses zero-based iteration indices, unlike HDF5 shell labels.
    params["keepall"] = ",".join(map(str, range(numerical["nmax"]))) if numerical.get("untruncated", False) else ""
    params["T"] = format(1e-5 * physical["gap"] / unit, ".17g")
    return "[param]\n" + "".join(f"{k}={v}\n" for k, v in params.items()) + f"[extra]\nB={physical['field'] / unit:.17g}\n"


def read_param(path) -> dict:
    """Validate a literal deck and return its two sections (case preserved)."""
    parser = ConfigParser(interpolation=None)
    parser.optionxform = str
    try:
        with Path(path).open() as stream:
            parser.read_file(stream)
    except ConfigError as error:
        raise ValueError(f"invalid param: {error}") from error
    if set(parser.sections()) != {"param", "extra"} or parser.defaults():
        raise ValueError("expected only [param] and [extra], without DEFAULT")
    params, extra = dict(parser["param"]), dict(parser["extra"])
    allowed = {*FIXED, *ENERGIES.values(), *NUMERICAL.values(), "bcsgap", "keepall", "T", "Nz"}
    if params.keys() - allowed or set(extra) != {"B"}:
        raise ValueError("unsupported param/extra keys (no spectra, polarization or alternate conventions)")
    for key, expected in FIXED.items():
        if params.get(key) != expected:
            raise ValueError(f"param {key} must be {expected!r}")
    try:
        physical = {key: float(params[deck]) for key, deck in ENERGIES.items()}
        physical.update(geometry="single", gap=physical["bandwidth"] * float(params["bcsgap"]), field=float(extra["B"]))
        numerical = {key: (int if key in ("nmax", "keep", "keepmin") else float)(params[deck])
                     for key, deck in NUMERICAL.items()}
        numerical["untruncated"] = bool(params["keepall"])
        render_param(physical, numerical)  # Validate ranges before expanding keepall.
        full = ",".join(map(str, range(numerical["nmax"])))
        if params["keepall"] not in ("", full):
            raise ValueError("keepall must be empty or list every zero-based iteration")
        if not math.isclose(float(params["T"]), 1e-5 * physical["gap"], rel_tol=1e-13, abs_tol=0):
            raise ValueError("T must equal 1e-5 times the gap in deck units")
        if "Nz" in params and int(params["Nz"]) < 1:
            raise ValueError("Nz must be a positive integer")
    except (KeyError, OverflowError) as error:
        raise ValueError(f"missing or invalid parameter: {error}") from error
    return {"param": params, "extra": extra}


def extract_result(directory) -> dict:
    """Read param/data/raw.h5, independent of controller receipts; never write.

    Return physical settings in canonical gap units, numerical settings, states,
    coverage and diagnostics. Retain all subgap states plus the first two sorted
    roots of every saved sector. Degenerate diagonal expectations are flagged,
    not asserted to be unique. No residual is inferred from an NRG spectrum.
    """
    import h5py
    import numpy as np

    if __package__:
        from .physics import finish_states, state_values
    else:
        from physics import finish_states, state_values

    directory = Path(directory)
    deck = read_param(directory / "param")
    params = deck["param"]
    gap = float(params["bandrescale"]) * float(params["bcsgap"])
    physical = {key: float(params[value]) / gap for key, value in ENERGIES.items()}
    physical.update(geometry="single", gap=1., field=float(deck["extra"]["B"]) / gap)
    numerical = {key: (int if key in ("nmax", "keep", "keepmin") else float)(params[value])
                 for key, value in NUMERICAL.items()}
    numerical["untruncated"] = bool(params["keepall"])
    iteration = numerical["nmax"]

    def real_array(value, label):
        array = np.asarray(value)
        if not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all() or np.any(np.abs(array.imag) > 1e-10):
            raise ValueError(f"{label}: expected finite real data")
        return array.real

    def dataset(group, path, ndim):
        if path not in group or not isinstance(group[path], h5py.Dataset):
            raise ValueError(f"raw.h5: missing dataset {group.name}/{path}")
        item = group[path]
        if (ndim == 0 and item.shape not in ((), (1,))) or (ndim != 0 and item.ndim != ndim):
            raise ValueError(f"raw.h5: invalid shape for {item.name}")
        return real_array(item[()], item.name)

    lines = [line.strip() for line in (directory / "data").read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    chain = {"source": "data: z/Z blocks", "n_bath_sites": iteration + 1, "n_bonds": iteration,
             "unused_last_bond_coefficient": True}
    for tag, names in (("z", ("xi", "zeta")), ("Z", ("delta", "kappa"))):
        if lines.count(tag) != 1:
            raise ValueError(f"data: expected one {tag} coefficient block")
        cursor = lines.index(tag) + 1
        for name in names:
            try:
                count = int(lines[cursor]) + 1
                if count != iteration + 1:
                    raise ValueError(f"data: wrong {name} coefficient count for Nmax")
                values = [float(lines[i].replace("D", "e").replace("d", "e")) / gap
                          for i in range(cursor + 1, cursor + 1 + count)]
            except IndexError as error:
                raise ValueError(f"data: incomplete {name} coefficient table") from error
            chain[name + "_over_gap"] = real_array(values, name).tolist()
            cursor += count + 1
    if (any(abs(d - 1) > 1e-12 for d in chain["delta_over_gap"])
            or any(abs(k) > 1e-12 for k in chain["kappa_over_gap"])):
        raise ValueError("bath vacuum requires uniform gap and no anomalous hopping")
    hopping = np.array(chain["xi_over_gap"][:-1])
    normal = np.diag(chain["zeta_over_gap"]) + np.diag(hopping, 1) + np.diag(hopping, -1)
    # NRG centers each bath onsite term at n-1, so there is no trace(normal).
    vacuum = float(-np.hypot(np.linalg.eigvalsh(normal), 1).sum())
    chain["impurity_hopping_over_gap"] = math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi)
    states, spectra = [], {}
    with h5py.File(directory / "raw.h5", "r") as h5:
        iterations = sorted(int(key) for key in h5 if key.isascii() and key.isdigit())
        if not iterations or iterations[-1] != iteration:
            raise ValueError("raw.h5: missing requested final Nmax iteration")
        raw_ground = float(dataset(h5, "stats/GS_energy", 0).item()) / gap
        offset = float(dataset(h5, f"stats/{iteration}/energy_offset", 0).item()) / gap
        steps = [int(key) for key in h5["stats"] if key.isascii() and key.isdigit()]
        if max(steps) != iteration or not math.isclose(offset, raw_ground, rel_tol=1e-11, abs_tol=1e-9):
            raise ValueError("final energy offset disagrees with GS_energy/final iteration")
        last = h5[str(iteration)]
        if not isinstance(last, h5py.Group) or "eigen" not in last or not isinstance(last["eigen"], h5py.Group):
            raise ValueError("raw.h5: missing final eigen spectrum")
        for sector, group in last["eigen"].items():
            if not re.fullmatch(r"-?(0|[1-9][0-9]*)", sector) or str(int(sector)) != sector:
                raise ValueError(f"invalid signed SPU1 SSZ sector: {sector}")
            if not isinstance(group, h5py.Group):
                raise ValueError("raw.h5: expected eigen sector group")
            relative = dataset(group, "absenergy_zero", 1) / gap
            absolute = dataset(group, "absenergy", 1) / gap
            if relative.shape != absolute.shape or not np.allclose(absolute, raw_ground + relative, rtol=1e-11, atol=1e-9):
                raise ValueError(f"inconsistent absolute/ground-relative energies in sector {sector}")
            count = len(relative)
            if not count:
                continue
            diagonal = {}
            for operator in ("n_d", "n_d_ud", "SZd", "SZ2d"):
                path = f"s/{operator}/{sector}/{sector}/matrix"
                matrix = dataset(last, path, 2)
                if matrix.shape != (count, count):
                    raise ValueError(f"{path}: matrix shape disagrees with eigenvalue count")
                if path + "-imag" in last:
                    imaginary = dataset(last, path + "-imag", 2)
                    if imaginary.shape != matrix.shape or np.any(np.abs(imaginary) > 1e-10):
                        raise ValueError(f"{path}: expected a real operator matrix")
                if not np.allclose(matrix, matrix.T, rtol=0, atol=1e-8):
                    raise ValueError(f"{path}: operator matrix is not Hermitian")
                diagonal[operator] = np.diag(matrix).copy()
            order = np.argsort(relative, kind="stable")
            sorted_relative = relative[order]
            twice_sz = int(sector) - 1
            spectra[sector] = {"twice_sz": twice_sz, "indices_in_energy_order": order.tolist(),
                               "ground_relative_energy_over_gap": sorted_relative.tolist()}
            for root, index in enumerate(order):
                values = state_values(*(float(diagonal[op][index]) for op in ("n_d", "n_d_ud", "SZd", "SZ2d")), twice_sz / 2)
                values["energy"] = float(absolute[index]) - vacuum - physical["detuning"]
                degenerate = ((root > 0 and sorted_relative[root] - sorted_relative[root - 1] <= 1e-9)
                              or (root + 1 < count and sorted_relative[root + 1] - sorted_relative[root] <= 1e-9))
                states.append({"twice_sz": twice_sz, "parity": twice_sz % 2, "root": root,
                               "values": values, "degenerate": bool(degenerate), "residual": None})
    if not spectra or abs(min(min(row["ground_relative_energy_over_gap"]) for row in spectra.values())) > 1e-9:
        raise ValueError("lowest all-sector ground-relative energy is not zero")
    ground = min(row["values"]["energy"] for row in states)
    states = [row for row in states if row["root"] < 2 or row["values"]["energy"] - ground < 1]
    states.sort(key=lambda row: (row["twice_sz"], row["root"]))
    coverage = finish_states(states)
    coverage.update(scope="all final-shell subgap roots and first two roots per saved sector",
                    same_sector_degeneracy_tolerance_over_gap=1e-9,
                    diagonal_observables_unique=not any(row["degenerate"] for row in states if row["subgap"]),
                    continuum_completeness_established=False)
    return {"physical": physical, "numerical": numerical, "states": states, "coverage": coverage,
            "diagnostics": {"param": deck, "deck_gap_scale": gap, "final_h5_iteration": iteration,
                            "wilson_chain": chain, "analytic_bath_vacuum_over_gap": vacuum,
                            "raw_total_energy_over_gap": raw_ground, "final_energy_offset_over_gap": offset,
                            "all_sector_spectra": spectra,
                            "energy_reference": "raw/gap - own centered bath vacuum/gap - detuning/gap",
                            "spin_convention": "SSZ=twice_sz+1; positive B raises spin up; scalar SZd/SZ2d",
                            "bath_observables": "fixed-M longitudinal identities, not transverse observables"}}
