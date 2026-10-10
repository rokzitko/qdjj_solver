"""Read-only two-lead SPU1 adapter; never launches NRG or Mathematica.

The paper deck needs calculation-local model.m and modeloperators.m snapshots.
COMPLEX preserves phase-dependent matrices (HDF5 real plus matrix-imag);
SCHUR selects the initializer's eigensystem method. Scalar operators are
recalculated in the first pass; no spectra are requested. The SDS representation is only a
WN=WS=0 control, not the production paper-basis model. Only extraction imports
NumPy/h5py; copying nrg.py and physics.py together supports standalone use.
"""

from configparser import ConfigParser, Error as ConfigError
import math
from pathlib import Path
import re


FIXED = {
    "symtype": "SPU1", "band": "flat", "Ninit": "0", "lastall": "true",
    "discretization": "Z", "tri": "old", "diagratio": "1", "diag_mode": "serial",
    "diagth": "1", "strategy": "kept", "finite": "true", "dm": "false",
    "fdm": "false", "broaden": "false", "savebins": "false", "calc0": "true",
    "h5raw": "true", "h5all": "false", "h5last": "true", "h5ops": "true",
    "h5vectors": "false", "h5U": "false", "removefiles": "true", "silent": "true",
}
LOCAL_OPS = ("n_d", "n_d_ud", "SZd", "SZ2d")
REPRESENTATIONS = {
    "paper": {"model": "model.m", "readdir": ".", "options": "COMPLEX SCHUR",
              "ops": " ".join((*LOCAL_OPS, "asqJ"))},
    "sds": {"model": "ONE", "variant": "SDSMAGFIELD", "options": "SCHUR",
            "ops": " ".join((*LOCAL_OPS, "hop0", "hop1"))},
}
NUMERICAL = {"lambda": "Lambda", "z": "z", "nmax": "Nmax", "keep": "keep",
             "keepmin": "keepmin", "keepenergy": "keepenergy"}
ENERGIES = {"bandwidth": "bandrescale", "u": "U", "gamma": "Gamma"}
PHYSICAL_KEYS = {*ENERGIES, "geometry", "gap", "nu", "field", "phi", "wn", "ws"}
LITERAL = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def render_param(physical, numerical, units="gap", representation="paper") -> str:
    """Render literal inputs, with total Gamma and dimensional WN/WS.

    All energies share the input physical unit. nu and phi are dimensionless;
    keepenergy is in Wilson-shell units. Lambda >= 4 is mandatory, including
    tiny untruncated controls. SDS needs U>0 to recover nu from its detuning.
    """
    if set(physical) != PHYSICAL_KEYS:
        raise ValueError("physical requires exactly geometry, gap, bandwidth, u, gamma, nu, field, phi, wn, ws")
    if not set(NUMERICAL) <= numerical.keys() or numerical.keys() - {*NUMERICAL, "untruncated"}:
        raise ValueError("unsupported or missing numerical settings")
    if units not in ("gap", "bandwidth") or physical["geometry"] != "junction":
        raise ValueError("units must be gap or bandwidth; geometry must be junction")
    if representation not in REPRESENTATIONS:
        raise ValueError("representation must be paper or sds")
    for label, values, keys in (("physical", physical, PHYSICAL_KEYS - {"geometry"}),
                                ("numerical", numerical, set(NUMERICAL))):
        for key in keys:
            value = values[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{label}.{key} must be finite and real")
    if min(physical["gap"], physical["bandwidth"]) <= 0 or min(physical["u"], physical["gamma"]) < 0:
        raise ValueError("gap/bandwidth must be positive; u/gamma must be nonnegative")
    if numerical["lambda"] < 4 or not 0 < numerical["z"] <= 1:
        raise ValueError("Lambda must be >= 4 and z in (0, 1]")
    for key, minimum in (("nmax", 1), ("keep", 2), ("keepmin", 0)):
        if type(numerical[key]) is not int or numerical[key] < minimum:
            raise ValueError(f"{key} must be an integer >= {minimum}")
    if numerical["nmax"] > 998 or numerical["keepmin"] > numerical["keep"]:
        raise ValueError("Nmax must be <= 998 and keepmin <= keep")
    if type(numerical.get("untruncated", False)) is not bool:
        raise ValueError("untruncated must be boolean")
    if representation == "sds" and (physical["wn"] != 0 or physical["ws"] != 0 or physical["u"] == 0):
        raise ValueError("SDS control requires WN=WS=0 and U>0 (to recover nu)")
    unit = physical[units]
    params = dict(FIXED, **REPRESENTATIONS[representation])
    params.update({deck: physical[key] / unit for key, deck in ENERGIES.items()})
    params.update({deck: numerical[key] for key, deck in NUMERICAL.items()})
    if representation == "paper":
        params["bcsgap"] = -physical["gap"] / physical["bandwidth"]
        extra = {"nu": physical["nu"], "Bz": physical["field"] / unit,
                 "WN": physical["wn"] / unit, "WS": physical["ws"] / unit, "phi": physical["phi"]}
    else:
        params.update(delta=physical["u"] * (1 - physical["nu"]) / unit,
                      bcsgap1=-physical["gap"] / physical["bandwidth"],
                      bcsgap2=physical["gap"] / physical["bandwidth"])
        extra = {"B": physical["field"] / unit, "phi": physical["phi"]}
    # keepall counts zero-based iterations, whereas HDF5 shell labels start at 1.
    params["keepall"] = ",".join(map(str, range(numerical["nmax"]))) if numerical.get("untruncated", False) else ""
    params["T"] = 1e-5 * (physical["gap"] / unit)
    for value in (*params.values(), *extra.values()):
        if isinstance(value, (int, float)) and not math.isfinite(value):
            raise ValueError("nonfinite parameter after unit conversion")
    return "[param]\n" + "".join(f"{k}={v if isinstance(v, str) else format(v, '.17g')}\n"
                                  for k, v in params.items()) + "[extra]\n" + "".join(
                                      f"{k}={float(v):.17g}\n" for k, v in extra.items())


def read_param(path) -> dict:
    """Validate the supported literal deck; include parsed internal metadata.

    physical is in deck units here, and in gap units in extract_result. Unknown
    fields, especially [extra] rules that could change physics, are rejected.
    """
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
    representation = "paper" if params.get("model") == "model.m" else "sds"
    fixed = dict(FIXED, **REPRESENTATIONS[representation])
    variable = {*ENERGIES.values(), *NUMERICAL.values(), "T"}
    variable |= {"bcsgap"} if representation == "paper" else {"delta", "bcsgap1", "bcsgap2"}
    extra_keys = {"nu", "Bz", "WN", "WS", "phi"} if representation == "paper" else {"B", "phi"}
    if params.keys() - {*fixed, *variable, "keepall", "Nz"} or set(extra) != extra_keys:
        raise ValueError("unsupported param/extra keys (no spectra, polarization or alternate conventions)")
    for key, expected in fixed.items():
        if params.get(key) != expected:
            raise ValueError(f"param {key} must be {expected!r}")
    try:
        for key in variable | ({"Nz"} & params.keys()):
            if not LITERAL.fullmatch(params[key]):
                raise ValueError(f"param {key} must be a literal real number")
        if any(not LITERAL.fullmatch(value) for value in extra.values()):
            raise ValueError("extra values must be literal real numbers")
        physical = {key: float(params[deck]) for key, deck in ENERGIES.items()}
        pairing = float(params["bcsgap" if representation == "paper" else "bcsgap1"])
        if not pairing < 0:
            raise ValueError("first pairing coefficient must be negative; the physical gap is positive")
        physical.update(geometry="junction", gap=abs(physical["bandwidth"] * pairing), phi=float(extra["phi"]))
        if representation == "paper":
            physical.update(nu=float(extra["nu"]), field=float(extra["Bz"]), wn=float(extra["WN"]), ws=float(extra["WS"]))
        else:
            if physical["u"] <= 0 or float(params["bcsgap2"]) != -pairing:
                raise ValueError("SDS requires U>0 and opposite equal pairing coefficients")
            physical.update(nu=1 - float(params["delta"]) / physical["u"], field=float(extra["B"]), wn=0., ws=0.)
        numerical = {key: (int if key in ("nmax", "keep", "keepmin") else float)(params[deck])
                     for key, deck in NUMERICAL.items()}
        numerical["untruncated"] = bool(params["keepall"])
        render_param(physical, numerical, representation=representation)
        if params["keepall"] not in ("", ",".join(map(str, range(numerical["nmax"])))):
            raise ValueError("keepall must be empty or list every zero-based iteration")
        if not math.isclose(float(params["T"]), 1e-5 * physical["gap"], rel_tol=1e-13, abs_tol=0):
            raise ValueError("T must equal 1e-5 times the gap in deck units")
        if "Nz" in params and int(params["Nz"]) < 1:
            raise ValueError("Nz must be a positive integer")
    except (KeyError, OverflowError) as error:
        raise ValueError(f"missing or invalid parameter: {error}") from error
    return {"param": params, "extra": extra, "representation": representation,
            "physical": physical, "numerical": numerical}


def extract_result(directory) -> dict:
    """Read final-shell states, retaining all subgap and first three sector roots.

    Subtract both disconnected centered bath vacua, then the representation's
    dot constant. current is <dH/dphi>/gap, in units 2 e Delta/hbar. Degenerate
    diagonal expectations are flagged as basis dependent, never hard-failed.
    """
    import h5py
    import numpy as np

    if __package__:
        from .physics import finish_states, state_values
    else:
        from physics import finish_states, state_values

    directory = Path(directory)
    deck = read_param(directory / "param")
    gap = deck["physical"]["gap"]
    physical = {key: value if key in ("geometry", "nu", "phi") else value / gap
                for key, value in deck["physical"].items()}
    numerical, representation = deck["numerical"], deck["representation"]
    iteration = numerical["nmax"]
    paper_offset = physical["u"] / 2 * (physical["nu"] - 1)**2
    dot_offset = paper_offset if representation == "paper" else physical["u"] * (1 - physical["nu"])

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
    chains = [{"source": "data: z/Z blocks", "channel": channel + 1,
               "n_bath_sites": iteration + 1, "n_bonds": iteration,
               "unused_last_bond_coefficient": True} for channel in range(2)]
    # Legacy tables are coefficient-major: xi1, xi2, zeta1, zeta2 (likewise Z).
    for tag, names in (("z", ("xi", "zeta")), ("Z", ("delta", "kappa"))):
        if lines.count(tag) != 1:
            raise ValueError(f"data: expected one {tag} coefficient block")
        cursor = lines.index(tag) + 1
        for name in names:
            for chain in chains:
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
        if cursor < len(lines) and LITERAL.fullmatch(lines[cursor]):
            raise ValueError(f"data: extra coefficient table in {tag}; expected two channels")
    vacuum = 0.
    for index, chain in enumerate(chains):
        pairing = -1 if representation == "paper" or index == 0 else 1
        if (any(abs(d - pairing) > 1e-12 for d in chain["delta_over_gap"])
                or any(abs(k) > 1e-12 for k in chain["kappa_over_gap"])):
            raise ValueError("bath vacuum requires the deck's uniform signed gaps and no anomalous hopping")
        hopping = np.array(chain["xi_over_gap"][:-1])
        normal = np.diag(chain["zeta_over_gap"]) + np.diag(hopping, 1) + np.diag(hopping, -1)
        # NRG centers onsite terms at n-1: no trace(normal) in either vacuum.
        vacuum += float(-np.hypot(np.linalg.eigvalsh(normal), 1).sum())
    gamma_total = math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi)
    states, spectra = [], {}
    with h5py.File(directory / "raw.h5", "r") as h5:
        if "params/channels" in h5 and dataset(h5, "params/channels", 0).item() != 2:
            raise ValueError("raw.h5: expected two channels")
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
            for operator in deck["param"]["ops"].split():
                path = f"s/{operator}/{sector}/{sector}/matrix"
                matrix = dataset(last, path, 2).astype(complex)
                if matrix.shape != (count, count):
                    raise ValueError(f"{path}: matrix shape disagrees with eigenvalue count")
                if path + "-imag" in last:
                    imaginary = dataset(last, path + "-imag", 2)
                    if imaginary.shape != matrix.shape:
                        raise ValueError(f"{path}: imaginary matrix shape mismatch")
                    matrix += 1j * imaginary
                if not np.allclose(matrix, matrix.conj().T, rtol=0, atol=1e-8):
                    raise ValueError(f"{path}: operator matrix is not Hermitian")
                diagonal[operator] = real_array(np.diag(matrix), path).copy()
            order = np.argsort(relative, kind="stable")
            sorted_relative = relative[order]
            twice_sz = int(sector) - 1
            spectra[sector] = {"twice_sz": twice_sz, "indices_in_energy_order": order.tolist(),
                               "ground_relative_energy_over_gap": sorted_relative.tolist()}
            for root, index in enumerate(order):
                values = state_values(*(float(diagonal[op][index]) for op in LOCAL_OPS), twice_sz / 2)
                values["energy"] = float(absolute[index]) - vacuum - dot_offset
                values["energy_paper"] = (float(absolute[index]) - vacuum if representation == "paper"
                                          else values["energy"] + paper_offset)
                values["current"] = (float(diagonal["asqJ"][index]) / gap if representation == "paper" else
                                     gamma_total / 4 * (-math.sin(physical["phi"] / 4) * float(diagonal["hop0"][index])
                                                        + math.cos(physical["phi"] / 4) * float(diagonal["hop1"][index])))
                degenerate = ((root > 0 and sorted_relative[root] - sorted_relative[root - 1] <= 1e-8)
                              or (root + 1 < count and sorted_relative[root + 1] - sorted_relative[root] <= 1e-8))
                states.append({"twice_sz": twice_sz, "parity": twice_sz % 2, "root": root,
                               "values": values, "degenerate": bool(degenerate), "residual": None})
    if not spectra or abs(min(min(row["ground_relative_energy_over_gap"]) for row in spectra.values())) > 1e-9:
        raise ValueError("lowest all-sector ground-relative energy is not zero")
    ground = min(row["values"]["energy"] for row in states)
    states = [row for row in states if row["root"] < 3 or row["values"]["energy"] - ground < 1]
    states.sort(key=lambda row: (row["twice_sz"], row["root"]))
    coverage = finish_states(states)
    coverage.update(scope="all final-shell subgap roots and first three roots per saved sector",
                    same_sector_degeneracy_tolerance_over_gap=1e-8,
                    diagonal_observables_unique=not any(row["degenerate"] for row in states if row["subgap"]),
                    continuum_completeness_established=False)
    return {"physical": physical, "numerical": numerical, "states": states, "coverage": coverage,
            "diagnostics": {"representation": representation, "deck_gap_scale": gap,
                            "final_h5_iteration": iteration, "wilson_chains": chains,
                            "analytic_bath_vacuum_over_gap": vacuum, "paper_offset_over_gap": paper_offset,
                            "raw_total_energy_over_gap": raw_ground, "final_energy_offset_over_gap": offset,
                            "all_sector_spectra": spectra,
                            "energy_reference": "raw/gap - two centered bath vacua/gap - " +
                                                ("U/(2 gap) (nu-1)^2" if representation == "paper" else "detuning/gap"),
                            "current_normalization": "<dH/dphi>/gap; I in units 2 e Delta/hbar",
                            "spin_convention": "SSZ=twice_sz+1; positive Bz raises spin up; scalar SZd/SZ2d",
                            "bath_observables": "fixed-M longitudinal identities, not transverse observables"}}
