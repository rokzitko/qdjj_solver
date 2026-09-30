"""Numerical definitions for the large-Gamma study (Delta = 1).

The quadratic references use physical-electron Nambu matrices and continuum
Green functions, independently of the many-body solvers' QP coordinates.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from qdjj_solver import Sector, reference_model
from qdjj_solver.common.io import fingerprint
from qdjj_solver.qp_solver import dimension

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "large_Gamma/input/study.json"
OBSERVABLES = ("impurity_charge", "double_occupancy_0", "impurity_spin_z", "phase_derivative")
VALUES = ("energy_even", "energy_odd", "signed_gap", "q_d", "P_0", "P_1", "P_2",
          "C_d_bath", "even_P_2", "current_even", "current_odd")
SECTORS = tuple(Sector(parity, parity, eta) for parity in (0, 1) for eta in (1, -1))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_config(config):
    expected = {"format_version", "model", "u_values", "gamma_values", "phases", "baths",
                "cutoffs", "chi_values", "qp_options", "dmrg_options", "convergence",
                "resources", "controls", "refinement"}
    if set(config) != expected or config["format_version"] != 1:
        raise ValueError("unsupported or incomplete large-Gamma configuration")
    model = config["model"]
    if model != dict(delta=1., bandwidth=100., detuning=0., field=0., rho_wn=0.,
                     rho_ws=0., symmetry=True, compress=True, bath_reference="isolated"):
        raise ValueError("this study requires the declared symmetric D=100, Delta=1 junction")
    for name in ("u_values", "gamma_values", "phases", "cutoffs", "chi_values"):
        values = config[name]
        if not values or not all(np.isfinite(v) for v in values) or len(set(values)) != len(values):
            raise ValueError(f"invalid {name}")
    if min(config["u_values"]) < 0 or min(config["gamma_values"]) <= 0:
        raise ValueError("nonnegative U and strictly positive Gamma required")
    if any(not 0 <= p < np.pi for p in config["phases"]):
        raise ValueError("phases must avoid the pi degeneracy")
    if any(q not in (2, 4, 6) for q in config["cutoffs"]):
        raise ValueError("the compared cutoffs are 2, 4, and 6")
    chis = config["chi_values"]
    if any(not isinstance(v, int) or v < 1 for v in chis) or sorted(chis) != chis:
        raise ValueError("bond limits must be strictly increasing positive integers")
    if config["convergence"]["stable_steps"] < 2:
        raise ValueError("at least two stable refinement steps required")
    for name in ("target", "safety_factor", "thresholds"):
        if np.any(np.asarray(config["convergence"][name]) <= 0):
            raise ValueError("convergence targets must be positive")
    for family in config["baths"]:
        if family["kind"] not in ("cosh-grid", "surrogate", "chain-expansion"):
            raise ValueError("unknown bath family")
        levels = family["levels"]
        if not levels or sorted(set(levels)) != levels or min(levels) < 2:
            raise ValueError("bath sizes must be increasing")
        if family["kind"] == "cosh-grid" and any(n % 2 for n in levels):
            raise ValueError("cosh grids need even signed-level counts")
    return config


def points(config):
    return [dict(u=float(u), gamma=float(g), phi=float(phi))
            for phi in config["phases"] for u in config["u_values"] for g in config["gamma_values"]]


def point_key(point):
    return fingerprint(point)[:20]


def make_hamiltonian(bath, point, config, *, compress=None):
    model = {k: v for k, v in config["model"].items() if k not in ("delta", "bandwidth")}
    if compress is not None:
        model["compress"] = compress
    return reference_model(bath, **model, **point)


def sector_key(sector):
    return f"p{sector.parity}_eta{sector.eta:+d}"


def values_from_sectors(sectors):
    """Select the lowest active-space state of each parity, comparing both etas.

    The doublet is followed even when excited. The free dark-channel threshold
    is stored separately; it must not replace this branch silently.
    """
    chosen = {p: min((s for s in sectors if s["sector"]["parity"] == p),
                     key=lambda s: s["energy"]) for p in (0, 1)}
    even, odd = chosen[0], chosen[1]
    obs = odd["observables"]
    charge, p2, qd = obs["impurity_charge"], obs["double_occupancy_0"], 2*obs["impurity_spin_z"]
    p1 = charge-2*p2
    gap = odd["energy"]-even["energy"]
    result = dict(energy_even=even["energy"], energy_odd=odd["energy"], signed_gap=gap,
                  q_d=qd, P_0=1-charge+p2, P_1=p1, P_2=p2, C_d_bath=.75*(qd-p1),
                  even_P_2=even["observables"]["double_occupancy_0"],
                  current_even=even["observables"]["phase_derivative"],
                  current_odd=obs["phase_derivative"])
    validate_values(result)
    return result


def validate_values(values):
    if set(values) != set(VALUES) or not all(np.isfinite(v) for v in values.values()):
        raise ValueError("missing/nonfinite physical measurements")
    probs = [values[f"P_{i}"] for i in range(3)]
    if (min(probs) < -2e-8 or max(probs) > 1+2e-8 or abs(sum(probs)-1) > 2e-8
            or abs(probs[0]-probs[2]) > 2e-8):
        raise ValueError("invalid particle-hole-symmetric charge probabilities")
    if not -values["P_1"]/3-2e-8 <= values["q_d"] <= values["P_1"]+2e-8:
        raise ValueError("state is incompatible with the assumed SU(2) doublet")


def gaussian_values(energy, excitation, projector, excited, response):
    n = len(projector)//2

    def local(p):
        up, down = p[0, 0].real, 1-p[n, n].real
        return dict(impurity_charge=float(up+down), impurity_spin_z=float((up-down)/2),
                    double_occupancy_0=float(up*down+abs(p[0, n])**2),
                    phase_derivative=float(np.trace(p @ response).real))

    return values_from_sectors([
        dict(sector=dict(parity=0), energy=energy, observables=local(projector)),
        dict(sector=dict(parity=1), energy=energy+excitation, observables=local(excited))])


def quadratic_finite(bath, point):
    """Exact U=0 even vacuum and spin-up odd quasiparticle, in electron space."""
    if point["u"] != 0 or point["gamma"] <= 0:
        raise ValueError("quadratic reference requires U=0 and Gamma>0")
    levels, gamma, phi = bath.levels, point["gamma"], point["phi"]
    n = 1+2*levels
    normal = np.diag(np.r_[0., bath.xi, bath.xi]).astype(complex)
    derivative = np.zeros_like(normal)
    pair = -np.diag(np.r_[0., np.full(2*levels, bath.delta)])
    for lead, velocity in enumerate((-.5, .5)):
        sl = slice(1+lead*levels, 1+(lead+1)*levels)
        hopping = np.sqrt(gamma/(2*np.pi*bath.rho)*bath.weights)*np.exp(.5j*velocity*phi)
        normal[sl, 0], normal[0, sl] = hopping, hopping.conj()
        derivative[sl, 0] = .5j*velocity*hopping
        derivative[0, sl] = derivative[sl, 0].conj()
    matrix = np.block([[normal, pair], [pair, -normal.conj()]])
    response = np.block([[derivative, np.zeros_like(pair)],
                         [np.zeros_like(pair), -derivative.conj()]])
    energies, vectors = np.linalg.eigh(matrix)
    occupied = vectors[:, :n]
    projector = occupied @ occupied.conj().T
    excited = projector+np.outer(vectors[:, n], vectors[:, n].conj())
    energy = float(energies[:n].sum()+2*bath.energies.sum())
    return gaussian_values(energy, float(energies[n]), projector, excited, response)


def quadratic_continuum(point, bandwidth=100.):
    """Independent finite-band integrals, including energy, Wick charge and spin.

    At half filling the positive bound-state spinor has equal electron/hole
    weights. Its dot residue Z is the inverse derivative of the pole equation.
    Occupying it changes the anomalous covariance F by -Z/2.
    """
    if point["u"] != 0 or point["gamma"] <= 0 or not abs(point["phi"]) < np.pi:
        raise ValueError("quadratic continuum reference outside its supported domain")
    gamma, phi = point["gamma"], point["phi"]
    cosine = np.cos(phi/2)

    def g(w):
        a = np.hypot(1., w)
        return 2/np.pi*np.arctan(bandwidth/a)/a

    def denominator(w):
        value = gamma*g(w)
        return (w*(1+value))**2+(cosine*value)**2

    errors = {}

    def integral(name, function):
        value, error = quad(function, 0., np.inf, epsabs=2e-11, epsrel=2e-12, limit=400)
        errors[name] = float(error)
        return value

    energy = -integral("energy", lambda w: np.log1p(
        2*gamma*g(w)+(gamma*g(w))**2+(cosine*gamma*g(w)/w)**2))/np.pi
    anomalous = gamma*cosine/np.pi*integral("anomalous", lambda w: g(w)/denominator(w))
    current = gamma**2*np.sin(phi)/(2*np.pi)*integral(
        "current", lambda w: g(w)**2/denominator(w))

    def retarded(e):
        a = np.sqrt((1-e)*(1+e))
        value = 2/np.pi*np.arctan(bandwidth/a)/a
        derivative = 2/np.pi*e*(bandwidth/(a*a*(a*a+bandwidth**2))
                                + np.arctan(bandwidth/a)/a**3)
        return value, derivative

    bound = brentq(lambda e: e+gamma*(e-cosine)*retarded(e)[0], 0., 1-1e-14, xtol=5e-15)
    value, derivative = retarded(bound)
    residue = 1/(1+gamma*value+gamma*(bound-cosine)*derivative)
    odd_current = current-.5*gamma*np.sin(phi/2)*value*residue
    p2_even = .25+anomalous**2
    p2 = p2_even-anomalous*residue
    values = dict(energy_even=float(energy), energy_odd=float(energy+bound), signed_gap=float(bound),
                  q_d=float(residue), P_0=float(p2), P_1=float(1-2*p2), P_2=float(p2),
                  C_d_bath=float(.75*(residue-1+2*p2)), even_P_2=float(p2_even),
                  current_even=float(current), current_odd=float(odd_current))
    validate_values(values)
    return dict(values=values, integration_errors=errors, bandwidth=bandwidth)


def refinement_changes(rows, stable_steps=2):
    """Max of the last independent refinement differences, per observable."""
    if len(rows) < stable_steps+1:
        return None
    recent = rows[-stable_steps-1:]
    return {key: max(abs(b["values"][key]-a["values"][key])
                     for a, b in zip(recent, recent[1:], strict=False)) for key in VALUES}


def empirical_resolution(rows, config, *, include_bond=True):
    changes = refinement_changes(rows, config["convergence"]["stable_steps"])
    if changes is None:
        return None
    recent = rows[-config["convergence"]["stable_steps"]-1:]
    if any(not row.get("accepted", False) for row in recent):
        return None
    safety = config["convergence"]["safety_factor"]
    return {k: max(config["convergence"]["target"], safety*(changes[k]+max(
        r.get("bond_changes", {}).get(k, 0.) if include_bond else 0. for r in recent))) for k in VALUES}


def resource_estimate(h, cutoff, sector, options):
    cutoff = h.bath_modes if cutoff is None else min(cutoff, h.bath_modes)
    d = dimension(h.nimp, h.spins, cutoff, sector, h.eta_labels)
    # Conservative complex Krylov allocation; sparse matrix is additional.
    ncv = min(d, options["ncv"])
    workspace = d*(8*((len(h.spins)+63)//64)+40+16*(2*ncv+9))
    return dict(dimension=d, workspace_bytes=workspace,
                allowed=d <= options["max_dimension"] and workspace < options["max_memory_gib"]*2**30)


def small_config():
    """A bounded, independent-reference validation fixture; never production data."""
    config = deepcopy(read_json(DEFAULT_INPUT))
    config.update(u_values=[0., 2.], gamma_values=[.4], phases=[1.8849555921538759],
                  chi_values=[8, 16, 32], baths=[dict(kind="cosh-grid", levels=[2, 4, 6])])
    config["dmrg_options"].update(seed_trials=2, max_sweeps=30)
    return config
