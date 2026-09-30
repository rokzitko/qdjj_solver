"""Static, single-reservoir NRG Ljubljana adapter; never launches a solver.

The installed nrginit/operators.m defines n_d, n_d_ud and the spin triplet
sigma_d. In initial.m, ireducsigma divides <S_z> by
CG(1/2,1/2;1,0 | 1/2,1/2) = 1/sqrt(3). Thus the doublet's up-member
moment is sigma_d_reduced/sqrt(3), NOT a Pauli expectation.

hamiltonian.m uses H_imp = delta*n + U/2*(n-1)^2, so its impurity energy
exceeds QD's by detuning. Its positive BCS pairing is related to QD's negative
pairing by a global gauge transformation; these static observables are invariant.
No bath vacuum or impurity constant is subtracted here.
"""

from configparser import ConfigParser
import math
from pathlib import Path


def _settings(physical, numerical, units):
    if units not in ("gap", "bandwidth"):
        raise ValueError("units must be 'gap' or 'bandwidth'")
    for key in ("gap", "bandwidth", "u", "gamma", "detuning", "field"):
        value = physical[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"physical.{key} must be a finite real number")
    if physical["gap"] <= 0 or physical["bandwidth"] <= 0:
        raise ValueError("gap and bandwidth must be positive")
    if physical["u"] < 0 or physical["gamma"] < 0:
        raise ValueError("u and total gamma must be nonnegative")
    if physical["geometry"] != "single" or physical["field"] != 0:
        raise ValueError("this SPSU2 adapter requires geometry='single' and field=0")
    for key in ("lambda", "z", "keepenergy"):
        value = numerical[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"numerical.{key} must be a finite real number")
    if numerical["lambda"] < 1.8:
        raise ValueError("Lambda < 1.8 is not authorized")
    if not 0 < numerical["z"] <= 1:
        raise ValueError("z must be in (0, 1]")
    for key, minimum in (("nmax", 1), ("keep", 2), ("keepmin", 0)):
        value = numerical[key]
        if type(value) is not int or value < minimum:
            raise ValueError(f"numerical.{key} must be an integer >= {minimum}")
    if numerical["nmax"] > 998:
        raise ValueError("nmax must not exceed NRG's limit of 998")
    if numerical["keepmin"] > numerical["keep"]:
        raise ValueError("keepmin must not exceed keep")
    if type(numerical.get("untruncated", False)) is not bool:
        raise ValueError("untruncated must be boolean")
    return float(physical[units])


def render_param(physical, numerical, units="gap", clean=False) -> str:
    """Render a literal [param] deck from physical energies in common units.

    Required numerical keys: lambda, z, nmax, keep, keepenergy, keepmin.
    There are no hidden truncation defaults. keepenergy is in Wilson-shell
    units, not physical units. Optional untruncated=True retains every state
    at every iteration; use only for small finite-chain controls.

    T is fixed at 1e-5*gap and scaled with the chosen energy unit. Nmax, not
    a temperature criterion, determines chain length. CLEAN removes the
    impurity, retaining exactly the same normal-band and pairing parameters.
    """
    unit = _settings(physical, numerical, units)
    if type(clean) is not bool:
        raise ValueError("clean must be boolean")
    gap = physical["gap"] / unit
    params = {
        "model": "CLEAN" if clean else "SIAM",
        "variant": "",
        "symtype": "SPSU2",
        "band": "flat",
        "bandrescale": physical["bandwidth"] / unit,
        # nrginit multiplies scdelta by bandrescale, just like xi and zeta.
        "bcsgap": physical["gap"] / physical["bandwidth"],
        "U": physical["u"] / unit,
        "Gamma": physical["gamma"] / unit,
        "delta": physical["detuning"] / unit,
        "Lambda": numerical["lambda"],
        "z": numerical["z"],
        "Ninit": 0,
        "Nmax": numerical["nmax"],
        "keep": numerical["keep"],
        "keepenergy": numerical["keepenergy"],
        "keepmin": numerical["keepmin"],
        # truncation.hpp tests ZERO-based iteration indices, not HDF5 labels.
        "keepall": ",".join(map(str, range(numerical["nmax"]))) if numerical.get("untruncated", False) else "",
        "lastall": True,
        "discretization": "Z",
        "tri": "old",
        "options": "SCHUR",
        "diagratio": 1,
        "diag_mode": "serial",
        "diagth": 1,
        "strategy": "kept",
        "ops": "" if clean else "n_d n_d_ud sigma_d",
        # Merely listing sigma_d in ops does not enable triplet recalculation.
        "spect": "" if clean else "sigma_d-sigma_d",
        "finite": not clean,
        "T": 1e-5 * gap,
        "dm": False,
        "fdm": False,
        "broaden": False,
        "savebins": False,
        "bins": 50,
        "calc0": True,
        "dumpannotated": 160,
        "dumpscaled": False,
        "dumpabs": False,
        "dumpEscale": gap,
        "dumpprecision": 16,
        "dumpgroups": False,
        "grouptol": 1e-10,
        "h5raw": True,
        "h5all": False,
        "h5last": True,
        "h5ops": True,
        "h5vectors": False,
        "h5U": False,
        "removefiles": True,
        "silent": True,
    }
    lines = ["[param]"]
    for key, value in params.items():
        if isinstance(value, bool):
            value = str(value).lower()
        elif isinstance(value, (float, int)):
            value = format(value, ".17g")
        lines.append(f"{key}={value}")
    return "\n".join(lines) + "\n"


def extract_result(workdir, physical, numerical, units="gap") -> dict:
    """Read raw.h5 and optional annotated.dat/data, returning JSON-safe values.

    All *_over_gap fields use the physical gap, independently of deck units.
    raw_total_energy fields retain the NRG deck units. branches contains the
    lowest singlet/doublet and their P0/P1/P2 and moment=<S_d^z> (up doublet).
    global_ground and all_sector_minima include every saved SPSU2 sector.
    Degenerate minima are flagged; individual expectations are not averaged.

    shell_convergence contains energy offsets/scales and annotated branch flow,
    not a claim of cutoff or continuum convergence. wilson_chain is None when
    data is absent; otherwise it contains the full generated coefficient tables
    in gap units. With Ninit=0, the finite chain uses sites 0..Nmax and bonds
    0..Nmax-1; the final generated xi/kappa entry is unused.

    Retain param for CLEAN identification (no impurity observables). Without
    param, SIAM is assumed and all impurity matrices are required. Missing or
    inconsistent required data raises rather than supplying zero observables.
    """
    unit = _settings(physical, numerical, units)
    import h5py  # Optional dependency: rendering and importing need no HDF5.

    workdir = Path(workdir)
    factor = unit / physical["gap"]
    model = "SIAM"
    if (workdir / "param").exists():
        parser = ConfigParser(interpolation=None)
        parser.optionxform = str
        parser.read(workdir / "param")
        params = parser["param"]
        model = params["model"]
        if model not in ("SIAM", "CLEAN") or params["symtype"] != "SPSU2" or params.get("variant", ""):
            raise ValueError("expected default SIAM or CLEAN with SPSU2 and empty variant")
        expected = ConfigParser(interpolation=None)
        expected.optionxform = str
        expected.read_string(render_param(physical, numerical, units, clean=model == "CLEAN"))
        for key in ("bandrescale", "bcsgap", "U", "Gamma", "delta", "Lambda", "z", "Ninit", "Nmax",
                    "keep", "keepenergy", "keepmin", "dumpEscale"):
            if not math.isclose(float(params[key]), float(expected["param"][key]), rel_tol=1e-13, abs_tol=0):
                raise ValueError(f"param {key} does not match the supplied manifest/units")
        for key in ("band", "discretization", "tri", "keepall", "lastall", "calc0",
                    "dumpscaled", "dumpabs", "dumpgroups"):
            if params[key] != expected["param"][key]:
                raise ValueError(f"param {key} is inconsistent with this adapter's deck")

    def real(value):
        value = complex(value)
        if not math.isfinite(value.real) or not math.isfinite(value.imag) or abs(value.imag) > 1e-10:
            raise ValueError("expected finite real NRG data")
        return float(value.real)

    def diagonal(last, path, index):
        value = complex(last[path][index, index])
        if path + "-imag" in last:
            value += 1j * real(last[path + "-imag"][index, index])
        return real(value)

    with h5py.File(workdir / "raw.h5", "r") as h5:
        iterations = sorted(int(key) for key in h5 if key.isdigit())
        if not iterations or iterations[-1] != numerical["nmax"]:
            raise ValueError("raw.h5 does not contain the requested final Nmax iteration")
        iteration = iterations[-1]
        last = h5[str(iteration)]
        raw_ground = real(h5["stats/GS_energy"][()])
        minima = {}
        for sector, group in last["eigen"].items():
            if not sector.isdigit() or int(sector) < 1:
                raise ValueError(f"not an SPSU2 spin-multiplicity sector: {sector}")
            relative = [real(e) * factor for e in group["absenergy_zero"][:]]
            absolute = [real(e) for e in group["absenergy"][:]]
            if len(relative) != len(absolute):
                raise ValueError(f"inconsistent eigenvalue counts in sector {sector}")
            if not relative:
                continue
            # absenergy_zero is shell-ground-relative, absenergy includes the
            # accumulated offset (c++/eigen.hpp). At the last shell it is GS.
            if any(not math.isclose(a * factor, raw_ground * factor + r, rel_tol=1e-11, abs_tol=1e-9)
                   for a, r in zip(absolute, relative, strict=True)):
                raise ValueError(f"inconsistent absolute/ground-relative energies in sector {sector}")
            index = min(range(len(relative)), key=relative.__getitem__)
            minima[sector] = {
                "spin_multiplicity": int(sector),
                "index": index,
                "ground_relative_energy_over_gap": relative[index],
                "raw_total_energy": absolute[index],
                "raw_total_energy_over_gap": absolute[index] * factor,
                "degenerate_minima": sum(abs(e - relative[index]) <= 1e-9 for e in relative),
                "computed_multiplets": len(relative),
            }
        if not {"1", "2"} <= minima.keys():
            raise ValueError("final spectrum must contain singlet and doublet sectors")
        ground = min(minima.values(), key=lambda row: row["ground_relative_energy_over_gap"])
        if abs(ground["ground_relative_energy_over_gap"]) > 1e-9:
            raise ValueError("lowest all-sector ground-relative energy is not zero")
        branches = {}
        for name, sector in (("singlet", "1"), ("doublet", "2")):
            branch = dict(minima[sector])
            if model == "SIAM":
                index = branch["index"]
                charge = diagonal(last, f"s/n_d/{sector}/{sector}/matrix", index)
                double = diagonal(last, f"s/n_d_ud/{sector}/{sector}/matrix", index)
                probabilities = {"P0": 1 - charge + double, "P1": charge - 2 * double, "P2": double}
                if any(p < -1e-8 or p > 1 + 1e-8 for p in probabilities.values()):
                    raise ValueError(f"unphysical impurity probabilities in {name}")
                reduced = diagonal(last, "t/sigma_d/2/2/matrix", index) if sector == "2" else 0.0
                moment = reduced / math.sqrt(3)
                if abs(moment) > probabilities["P1"] / 2 + 1e-8:
                    raise ValueError("impurity moment exceeds singly occupied probability / 2")
                branch.update(n_d=charge, n_d_ud=double, **probabilities, moment=moment)
                if sector == "2":
                    branch["sigma_d_reduced"] = reduced
            branches[name] = branch
        stats = []
        for step in sorted(int(key) for key in h5["stats"] if key.isdigit()):
            group = h5[f"stats/{step}"]
            row = {"iteration": step}
            for key in ("energyscale", "scale", "Teff", "abs_Egs", "energy_offset"):
                row[key + "_over_gap"] = real(group[key][()]) * factor
            stats.append(row)
        if not stats or stats[-1]["iteration"] != iteration:
            raise ValueError("missing final-shell statistics")
        if not math.isclose(stats[-1]["energy_offset_over_gap"], raw_ground * factor, rel_tol=1e-11, abs_tol=1e-9):
            raise ValueError("final energy offset disagrees with GS_energy")

    flows = []
    annotated = workdir / "annotated.dat"
    if annotated.exists():
        blocks = annotated.read_text().strip().split("\n\n")
        blocks = [block for block in blocks if block.strip()]
        if len(blocks) != iteration + 1:
            raise ValueError("annotated.dat must contain the initial cluster and every NRG shell")
        for step, block in enumerate(blocks):
            levels = {}
            for line in block.splitlines():
                if not line.strip():
                    continue
                energy, sector = line.split()
                if not sector.isdigit() or int(sector) < 1:
                    raise ValueError("annotated.dat is not ungrouped SPSU2 output")
                # dumpEscale already puts annotated energies in gap units.
                energy = real(float(energy))
                levels[sector] = min(energy, levels.get(sector, math.inf))
            flows.append({"iteration": step, "sector_minima_over_gap": levels,
                          "doublet_minus_singlet_over_gap": levels["2"] - levels["1"] if {"1", "2"} <= levels.keys() else None})
        for sector, energy in flows[-1]["sector_minima_over_gap"].items():
            if sector not in minima or abs(energy - minima[sector]["ground_relative_energy_over_gap"]) > 1e-9:
                raise ValueError("final annotated.dat energies disagree with raw.h5")

    chain = None
    data = workdir / "data"
    if data.exists():
        lines = [line.strip() for line in data.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        chain = {"source": "data: z/Z blocks", "n_bath_sites": iteration + 1,
                 "n_bonds": iteration, "unused_last_bond_coefficient": True}
        for tag, names in (("z", ("xi", "zeta")), ("Z", ("delta", "kappa"))):
            if lines.count(tag) != 1:
                raise ValueError(f"expected one {tag} Wilson-coefficient block in data")
            cursor = lines.index(tag) + 1
            for name in names:
                count = int(lines[cursor]) + 1
                if count != iteration + 1:
                    raise ValueError(f"wrong {name} coefficient count for Nmax")
                values = lines[cursor + 1:cursor + 1 + count]
                if len(values) != count:
                    raise ValueError(f"incomplete {name} coefficient table")
                chain[name + "_over_gap"] = [real(float(v.replace("D", "e").replace("d", "e"))) * factor for v in values]
                cursor += count + 1
        # Flat band: theta=2 and gammaA=1 for discretization=Z (nrginit).
        # This contact is analytic, unlike the generated chain coefficients.
        chain["impurity_hopping_over_gap"] = 0.0 if model == "CLEAN" else math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi) / physical["gap"]
        chain["impurity_hopping_source"] = "sqrt(2*D*totalGamma/pi), flat band, discretization Z"

    changes = {}
    for lag, name in ((1, "last_step_gap_change_over_gap"), (2, "same_parity_gap_change_over_gap")):
        changes[name] = None
        if len(flows) > lag:
            now, before = (flows[i]["doublet_minus_singlet_over_gap"] for i in (-1, -1 - lag))
            if now is not None and before is not None:
                changes[name] = now - before
    return {
        "backend": "nrg", "model": model, "units": units, "unit_value": unit,
        "physical": dict(physical), "numerical": dict(numerical),
        "final_h5_iteration": iteration,
        "raw_total_energy": raw_ground, "raw_total_energy_over_gap": raw_ground * factor,
        "branches": branches,
        "doublet_minus_singlet_over_gap": minima["2"]["ground_relative_energy_over_gap"] - minima["1"]["ground_relative_energy_over_gap"],
        "global_ground": dict(ground), "all_sector_minima": minima,
        "ground_spin_multiplicities": sorted(int(ss) for ss, row in minima.items() if abs(row["ground_relative_energy_over_gap"]) <= 1e-9),
        "shell_convergence": {"statistics": stats, "branch_flow": flows, **changes},
        "wilson_chain": chain,
        "energy_reference": "NRG total, no bath subtraction; H_imp=detuning*n+U/2*(n-1)^2",
        "moment_definition": "<S_d^z> in the S_z=+1/2 doublet; sigma_d_reduced/sqrt(3)",
    }
