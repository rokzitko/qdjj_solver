"""Standalone ASQ junction conventions and longitudinal state records (stdlib only)."""

import math
from numbers import Number


PHYSICAL = dict(geometry="junction", gap=1., bandwidth=100., u=2., gamma=.4,
                nu=.75, field=.03, phi=3*math.pi/5, wn=0., ws=0.)
PRIMARY = {0: "even", -1: "odd_down", 1: "odd_up"}
OBSERVABLES = ("energy", "P0", "P1", "P2", "P_up", "P_down", "m_dot", "m_bath",
               "dot_z2", "bath_z2", "Czz", "Czz_connected", "energy_paper", "current")
TOLERANCE = 1e-8


def _real(value, name):
    if not isinstance(value, Number) or isinstance(value, bool):
        raise ValueError(f"{name} must be finite and real")
    z = complex(value)
    if not math.isfinite(z.real) or not math.isfinite(z.imag) or abs(z.imag) > TOLERANCE:
        raise ValueError(f"{name} must be finite and real")
    return float(z.real)


def detuning(physical):
    """Centered gate energy U(1-nu), including the well-defined U=0 limit."""
    return physical["u"]*(1-physical["nu"])


def paper_offset(physical):
    """H_paper = H_centered + U(nu-1)^2/2, in dimensional energy units."""
    return physical["u"]/2*(physical["nu"]-1)*(physical["nu"]-1)


def validate_physical(physical):
    """Return a validated copy; nu is the only accepted gate specification.

    Gamma is total hybridization; W_N/W_S, unlike nu/phi, are energies.
    No bound on finite nu is imposed and the caller's input is never modified.
    """
    if set(physical) != set(PHYSICAL):
        raise ValueError("physical requires exactly PHYSICAL keys; nu is the authoritative gate")
    if physical["geometry"] != "junction":
        raise ValueError("test4 requires physical L/R junction geometry")
    checked = {name: _real(value, name) if name != "geometry" else value
               for name, value in physical.items()}
    if checked["gap"] <= 0 or checked["bandwidth"] <= 0 or checked["gamma"] < 0 or checked["u"] < 0:
        raise ValueError("gap/bandwidth must be positive and u/gamma nonnegative")
    _real(detuning(checked), "derived detuning")
    _real(paper_offset(checked), "paper energy offset")
    return checked


def state_values(charge, double, md, zd2, M, *, m_bath=None, Czz=None, bath_z2=None):
    """Return probabilities and longitudinal moments, without energies/current.

    P_up/P_down are singly occupied dot probabilities, not electron densities.
    Direct measurements of the total L+R bath spin are retained. Fixed-M
    identities use only spin U(1), never an assumption about total S^2.
    """
    charge, double, md, zd2, M = (
        _real(value, name) for name, value in
        zip(("charge", "double", "md", "zd2", "M"), (charge, double, md, zd2, M), strict=True)
    )
    p0, p1, p2 = 1-charge+double, charge-2*double, double
    probabilities = dict(P0=p0, P1=p1, P2=p2, P_up=p1/2+md, P_down=p1/2-md)
    if any(not -TOLERANCE <= p <= 1+TOLERANCE for p in probabilities.values()):
        raise ValueError("dot probabilities outside [0, 1]")
    if abs(p0+p1+p2-1) > TOLERANCE or abs(zd2-p1/4) > TOLERANCE:
        raise ValueError("dot normalization or dot_z2=P1/4 identity failed")
    direct = {"m_bath": m_bath, "Czz": Czz, "bath_z2": bath_z2}
    expected = {"m_bath": M-md, "Czz": M*md-zd2, "bath_z2": M*M-2*M*md+zd2}
    for name, value in direct.items():
        direct[name] = expected[name] if value is None else _real(value, name)
        if abs(direct[name]-expected[name]) > TOLERANCE:
            raise ValueError(f"fixed-M identity failed for {name}")
    if (not -TOLERANCE <= zd2 <= .25+TOLERANCE or zd2-md*md < -TOLERANCE
            or direct["bath_z2"]-direct["m_bath"]**2 < -TOLERANCE):
        raise ValueError("longitudinal moment bounds failed")
    connected = direct["Czz"]-md*direct["m_bath"]
    if abs(connected-(md*md-zd2)) > TOLERANCE:
        raise ValueError("connected fixed-M identity failed")
    values = dict(probabilities, m_dot=md, dot_z2=zd2, **direct, Czz_connected=connected)
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("derived observables must be finite")
    return values


def finish_states(states):
    """Validate/mutate all rows and report selected-state coverage, not completeness.

    Energies/excitations are in gap units; current is I/(2e Delta/hbar), real.
    Cross-M ties are symmetry resolved, not same-sector eigenvector ambiguities.
    Three subgap states are only a hypothesis; extra even/odd roots are retained.
    """
    if not states:
        raise ValueError("at least one state is required")
    seen, offset = set(), None
    for row in states:
        spin, parity, root = row["twice_sz"], row["parity"], row["root"]
        if (type(spin) is not int or type(parity) is not int or parity != spin % 2
                or type(root) is not int or root < 0):
            raise ValueError("invalid state sector or root")
        if (spin, root) in seen:
            raise ValueError("duplicate sector/root")
        seen.add((spin, root))
        if type(row["degenerate"]) is not bool:
            raise ValueError("degenerate must be boolean")
        values = row["values"]
        if set(values) != set(OBSERVABLES):
            raise ValueError("state values must contain all and only OBSERVABLES")
        values.update({name: _real(value, name) for name, value in values.items()})
        difference = values["energy_paper"]-values["energy"]
        if offset is not None and abs(difference-offset) > TOLERANCE:
            raise ValueError("paper energy offset must be constant across states")
        offset = difference if offset is None else offset
        checked = state_values(values["P1"]+2*values["P2"], values["P2"], values["m_dot"],
                               values["dot_z2"], spin/2, m_bath=values["m_bath"],
                               Czz=values["Czz"], bath_z2=values["bath_z2"])
        if any(abs(values[name]-value) > TOLERANCE for name, value in checked.items()):
            raise ValueError("inconsistent state probabilities or connected correlation")
        if row["residual"] is not None:
            row["residual"] = _real(row["residual"], "residual")
            if row["residual"] < 0:
                raise ValueError("residual must be nonnegative")
        row["id"] = PRIMARY[spin] if root == 0 and spin in PRIMARY else f"sz{spin:+d}-r{root}"
    ground = min(row["values"]["energy"] for row in states)
    for row in states:
        energy = row["values"]["energy"]
        row["excitation"] = energy-ground
        row["subgap"] = row["excitation"] < 1.
        row["degenerate"] |= any(other is not row and other["twice_sz"] == row["twice_sz"] and
                                 abs(other["values"]["energy"]-energy) <= TOLERANCE
                                 for other in states)
    subgap = [row["id"] for row in states if row["subgap"]]
    return dict(subgap_ids=subgap,
                degenerate_subgap_ids=[row["id"] for row in states if row["subgap"] and row["degenerate"]],
                ground_ids=[row["id"] for row in states if row["excitation"] <= TOLERANCE],
                sectors=sorted({row["twice_sz"] for row in states}),
                expected_three_subgap_states=set(subgap) == set(PRIMARY.values()),
                expected_three_subgap_states_is_hypothesis=True,
                paper_energy_offset=offset, current_units="I/(2e*Delta/hbar)", current_is_real=True,
                ground_energy=ground, scope="selected sectors and returned roots only")
