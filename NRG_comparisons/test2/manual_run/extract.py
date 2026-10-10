"""Standalone, read-only extraction for the fixed test2 SIAM/SPSU2 calculation.

Energies use Delta=1, D=100, U=2, total Gamma=0.4 and detuning=0.5.
Only extraction needs NumPy and h5py. No solver, repository imports, completion
markers or freshness receipts are used here; the controller owns run identity.
"""

from configparser import ConfigParser, Error as ConfigError
import math
from pathlib import Path


def read_param(path) -> dict[str, str]:
    """Validate a literal [param] deck, preserving all keys, case and values.

    Lambda, Nmax, keep, keepmin and keepenergy are editable numerical settings.
    z is in (0, 1]; Nz, when present, is a positive integer. Cross-directory
    agreement of every parameter except z must be checked by the caller.
    """
    parser = ConfigParser(interpolation=None)
    parser.optionxform = str
    try:
        with Path(path).open() as stream:
            parser.read_file(stream)
    except ConfigError as error:
        raise ValueError(f"invalid param file {path}: {error}") from error
    if parser.sections() != ["param"] or parser.defaults():
        raise ValueError("param requires one [param] section and no [DEFAULT] values")
    params = dict(parser["param"])
    fixed = {
        "model": "SIAM", "variant": "", "symtype": "SPSU2", "band": "flat",
        "lastall": "true", "keepall": "", "discretization": "Z", "tri": "old",
        "options": "SCHUR", "strategy": "kept", "ops": "n_d n_d_ud sigma_d",
        "spect": "sigma_d-sigma_d", "finite": "true", "dm": "false", "fdm": "false",
        "calc0": "true", "dumpscaled": "false", "dumpabs": "false", "dumpgroups": "false",
        "h5raw": "true", "h5all": "false", "h5last": "true", "h5ops": "true",
        "h5vectors": "false", "h5U": "false",
        "diag_mode": "serial", "broaden": "false", "savebins": "false",
        "removefiles": "true", "silent": "true",
    }
    fixed_numbers = {
        "bandrescale": 100., "bcsgap": .01, "U": 2., "Gamma": .4, "delta": .5,
        "Ninit": 0., "diagratio": 1., "dumpEscale": 1., "T": 1e-5,
        "diagth": 1., "bins": 50., "dumpannotated": 160., "dumpprecision": 16., "grouptol": 1e-10,
    }
    editable = {"Lambda", "Nmax", "keep", "keepmin", "keepenergy", "z", "Nz", "B"}
    unsupported = params.keys() - fixed.keys() - fixed_numbers.keys() - editable
    if unsupported:
        raise ValueError(f"unsupported param keys for fixed test2 calculation: {sorted(unsupported)}")
    for key, expected in fixed.items():
        if params.get(key) != expected:
            raise ValueError(f"param {key} must be {expected!r} for this adapter")

    def number(key, *, integer=False):
        try:
            value = int(params[key]) if integer else float(params[key])
        except (KeyError, ValueError) as error:
            kind = "integer" if integer else "finite number"
            raise ValueError(f"param {key} must be a {kind}") from error
        if not integer and not math.isfinite(value):
            raise ValueError(f"param {key} must be finite")
        return value

    for key, expected in fixed_numbers.items():
        if not math.isclose(number(key), expected, rel_tol=1e-13, abs_tol=0):
            raise ValueError(f"param {key} must be {expected:g} for fixed test2 physics/conventions")
    if "B" in params and number("B") != 0:
        raise ValueError("param B must be zero for the zero-field SPSU2 calculation")
    if number("Lambda") < 1.8:
        raise ValueError("param Lambda must be >= 1.8")
    if not 0 < number("z") <= 1:
        raise ValueError("param z must be in (0, 1]")
    if not 1 <= number("Nmax", integer=True) <= 998:
        raise ValueError("param Nmax must be an integer in 1..998")
    keep = number("keep", integer=True)
    if keep < 2:
        raise ValueError("param keep must be >= 2")
    if not 0 <= number("keepmin", integer=True) <= keep:
        raise ValueError("param keepmin must be between zero and keep")
    number("keepenergy")
    if "Nz" in params and number("Nz", integer=True) < 1:
        raise ValueError("param Nz must be a positive integer")
    return params


def extract_result(workdir) -> dict:
    """Return {values, diagnostics} from param, data and raw.h5, without writes.

    values has exactly ten entries: signed_gap, both branches' energy/P0/P1/P2,
    and doublet.moment. Energies subtract this chain's centered BCS vacuum and
    the NRG impurity constant 0.5. The moment is <S_d^z> of the up doublet,
    not a Pauli expectation. Repeated lowest same-sector multiplets are rejected
    because their individual expectations can depend on the eigenvector basis.
    Validation tolerances match the shared test1 adapter (not accuracy claims).
    """
    workdir = Path(workdir)
    params = read_param(workdir / "param")
    iteration = int(params["Nmax"])
    for name in ("data", "raw.h5"):
        if not (workdir / name).is_file():
            raise ValueError(f"missing required NRG output: {workdir / name}")

    import h5py
    import numpy as np

    def real(value, label):
        try:
            value = complex(value)
        except (ValueError, TypeError) as error:
            raise ValueError(f"{label}: expected finite real NRG data") from error
        if not math.isfinite(value.real) or not math.isfinite(value.imag) or abs(value.imag) > 1e-10:
            raise ValueError(f"{label}: expected finite real NRG data")
        return float(value.real)

    def dataset(group, path, ndim):
        if not isinstance(group, h5py.Group) or path not in group or not isinstance(group[path], h5py.Dataset):
            raise ValueError(f"raw.h5: missing dataset {group.name}/{path}")
        result = group[path]
        if ndim == 0:
            # NRG also stores scalar statistics as one-element vectors.
            if result.shape not in ((), (1,)):
                raise ValueError(f"raw.h5: {result.name} must be a scalar or one-element vector")
        elif result.ndim != ndim:
            raise ValueError(f"raw.h5: {result.name} must have {ndim} dimensions")
        return result

    def diagonal(group, path, index, count):
        matrix = dataset(group, path, 2)
        if matrix.shape != (count, count):
            raise ValueError(f"raw.h5: {matrix.name} matrix shape disagrees with eigenvalue count")
        value = complex(matrix[index, index])
        if path + "-imag" in group:
            imaginary = dataset(group, path + "-imag", 2)
            if imaginary.shape != matrix.shape:
                raise ValueError(f"raw.h5: {imaginary.name} matrix shape disagrees with real matrix")
            value += 1j * real(imaginary[index, index], imaginary.name)
        return real(value, matrix.name)

    with h5py.File(workdir / "raw.h5", "r") as h5:
        iterations = sorted(int(key) for key in h5 if key.isdigit())
        if not iterations or iterations[-1] != iteration or str(iteration) not in h5:
            raise ValueError("raw.h5 does not contain the requested final Nmax iteration")
        last = h5[str(iteration)]
        raw_ground = real(dataset(h5, "stats/GS_energy", 0)[()].item(), "stats/GS_energy")
        if not isinstance(last, h5py.Group) or "eigen" not in last or not isinstance(last["eigen"], h5py.Group):
            raise ValueError("raw.h5: missing final eigen spectrum")
        minima = {}
        for sector, group in last["eigen"].items():
            if not sector.isascii() or not sector.isdigit() or int(sector) < 1:
                raise ValueError(f"not an SPSU2 spin-multiplicity sector: {sector}")
            relative = [real(e, f"sector {sector} absenergy_zero")
                        for e in dataset(group, "absenergy_zero", 1)[:]]
            absolute = [real(e, f"sector {sector} absenergy")
                        for e in dataset(group, "absenergy", 1)[:]]
            if len(relative) != len(absolute):
                raise ValueError(f"inconsistent eigenvalue counts in sector {sector}")
            if not relative:
                continue
            if any(not math.isclose(a, raw_ground + r, rel_tol=1e-11, abs_tol=1e-9)
                   for a, r in zip(absolute, relative, strict=True)):
                raise ValueError(f"inconsistent absolute/ground-relative energies in sector {sector}")
            index = min(range(len(relative)), key=relative.__getitem__)
            degeneracy = sum(abs(e - relative[index]) <= 1e-9 for e in relative)
            if sector in {"1", "2"} and degeneracy > 1:
                raise ValueError(f"degenerate lowest multiplets in sector {sector}; expectations are basis dependent")
            minima[sector] = {
                "spin_multiplicity": int(sector), "index": index,
                "ground_relative_energy_over_gap": relative[index],
                "raw_total_energy_over_gap": absolute[index],
                "computed_multiplets": len(relative),
                "degenerate_minima": degeneracy,
            }
        if not {"1", "2"} <= minima.keys():
            raise ValueError("final spectrum must contain singlet and doublet sectors")
        ground = min(minima.values(), key=lambda row: row["ground_relative_energy_over_gap"])
        if abs(ground["ground_relative_energy_over_gap"]) > 1e-9:
            raise ValueError("lowest all-sector ground-relative energy is not zero")
        steps = sorted(int(key) for key in h5["stats"] if key.isdigit())
        if not steps or steps[-1] != iteration:
            raise ValueError("missing final-shell statistics")
        offset = real(dataset(h5, f"stats/{iteration}/energy_offset", 0)[()].item(), "final energy_offset")
        if not math.isclose(offset, raw_ground, rel_tol=1e-11, abs_tol=1e-9):
            raise ValueError("final energy offset disagrees with GS_energy")
        branches = {}
        for name, sector in (("singlet", "1"), ("doublet", "2")):
            branch = dict(minima[sector])
            index, count = branch["index"], branch["computed_multiplets"]
            charge = diagonal(last, f"s/n_d/{sector}/{sector}/matrix", index, count)
            double = diagonal(last, f"s/n_d_ud/{sector}/{sector}/matrix", index, count)
            probabilities = {"P0": 1 - charge + double, "P1": charge - 2 * double, "P2": double}
            if any(p < -1e-8 or p > 1 + 1e-8 for p in probabilities.values()):
                raise ValueError(f"unphysical impurity probabilities in {name}")
            branch.update(n_d=charge, n_d_ud=double, **probabilities)
            if sector == "2":
                reduced = diagonal(last, "t/sigma_d/2/2/matrix", index, count)
                moment = reduced / math.sqrt(3)
                if abs(moment) > probabilities["P1"] / 2 + 1e-8:
                    raise ValueError("impurity moment exceeds singly occupied probability / 2")
                branch.update(sigma_d_reduced=reduced, moment=moment)
            branches[name] = branch

    lines = [line.strip() for line in (workdir / "data").read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    chain = {"source": "data: z/Z blocks", "n_bath_sites": iteration + 1,
             "n_bonds": iteration, "unused_last_bond_coefficient": True}
    for tag, names in (("z", ("xi", "zeta")), ("Z", ("delta", "kappa"))):
        if lines.count(tag) != 1:
            raise ValueError(f"expected one {tag} Wilson-coefficient block in data")
        cursor = lines.index(tag) + 1
        for name in names:
            try:
                count = int(lines[cursor]) + 1
            except (IndexError, ValueError) as error:
                raise ValueError(f"data: missing or malformed {name} coefficient count") from error
            if count != iteration + 1:
                raise ValueError(f"data: wrong {name} coefficient count for Nmax")
            table = lines[cursor + 1:cursor + 1 + count]
            if len(table) != count:
                raise ValueError(f"data: incomplete {name} coefficient table")
            chain[name + "_over_gap"] = [
                real(v.replace("D", "e").replace("d", "e"), f"data: {name} coefficient") for v in table
            ]
            cursor += count + 1
    if (any(abs(d - 1) > 1e-12 for d in chain["delta_over_gap"])
            or any(abs(k) > 1e-12 for k in chain["kappa_over_gap"])):
        raise ValueError("vacuum subtraction requires uniform unit gap and no anomalous hopping")
    hopping = np.array(chain["xi_over_gap"][:-1])
    normal = np.diag(chain["zeta_over_gap"]) + np.diag(hopping, 1) + np.diag(hopping, -1)
    # NRG onsite energies are centered: no trace(normal) term belongs here.
    vacuum = real(-np.hypot(np.linalg.eigvalsh(normal), 1).sum(), "centered bath vacuum")
    chain["impurity_hopping_over_gap"] = math.sqrt(80 / math.pi)
    values = {}
    for name, branch in branches.items():
        values[f"{name}.energy"] = real(branch["raw_total_energy_over_gap"] - vacuum - .5,
                                        f"{name} corrected energy")
        values.update({f"{name}.{key}": branch[key] for key in ("P0", "P1", "P2")})
    values["doublet.moment"] = branches["doublet"]["moment"]
    gap = branches["doublet"]["ground_relative_energy_over_gap"] - branches["singlet"]["ground_relative_energy_over_gap"]
    if not math.isclose(values["doublet.energy"] - values["singlet.energy"], gap, rel_tol=1e-11, abs_tol=1e-9):
        raise ValueError("corrected branch difference disagrees with ground-relative signed gap")
    values["signed_gap"] = gap

    annotated = workdir / "annotated.dat"
    if annotated.exists():
        blocks = [block for block in annotated.read_text().strip().split("\n\n") if block.strip()]
        if len(blocks) != iteration + 1:
            raise ValueError("annotated.dat must contain the initial cluster and every NRG shell")
        levels = {}
        for line in blocks[-1].splitlines():
            fields = line.split()
            if len(fields) != 2 or fields[1] not in minima:
                raise ValueError("annotated.dat: expected ungrouped saved SPSU2 sectors")
            energy, sector = fields
            levels[sector] = min(real(energy, "annotated.dat energy"), levels.get(sector, math.inf))
        if any(abs(e - minima[s]["ground_relative_energy_over_gap"]) > 1e-9 for s, e in levels.items()):
            raise ValueError("final annotated.dat energies disagree with raw.h5")

    return {"values": values, "diagnostics": {
        "param": params, "final_h5_iteration": iteration,
        "raw_total_energy_over_gap": raw_ground, "final_energy_offset_over_gap": offset,
        "branches": branches, "all_sector_minima": minima, "global_ground": dict(ground),
        "ground_spin_multiplicities": sorted(int(s) for s, row in minima.items()
                                             if abs(row["ground_relative_energy_over_gap"]) <= 1e-9),
        "wilson_chain": chain, "analytic_bath_vacuum_over_gap": vacuum,
        "energy_reference": "raw NRG - own centered BCS bath vacuum - detuning (0.5)",
        "moment_definition": "<S_d^z> in the up doublet; sigma_d_reduced/sqrt(3)",
    }}
