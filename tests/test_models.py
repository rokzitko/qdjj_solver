"""Small physical-electron oracles for impurity and reservoir coordinates."""

from itertools import combinations, product

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.sparse import csr_matrix

from qdjj_solver.common.algebra import FermionOperator, annihilate, create, number
from qdjj_solver.common.baths import DiscreteBath
from qdjj_solver.common.models import Impurity, make_model, reference_model
from qdjj_solver.common.problem import Sector
from electron_oracle import electron_hamiltonian, electron_operators, evaluate, sector_indices


def linear_coordinates(model):
    records = model.metadata["coordinates"]["electron_annihilators"]
    u = np.zeros((len(records), len(model.spins)), dtype=complex)
    v = np.zeros_like(u)
    for row, record in enumerate(records):
        for term in record:
            (index,) = term["operators"]
            (u if index < 0 else v)[row, abs(index) - 1] = complex(term["real"], term["imag"])
    return u, v


def electron_coordinate_unitary(model, cs):
    """Build the many-body change of basis from a vacuum and CAR, independently.

    The inverse canonical map is a_j=sum_i(U_ij* c_i+V_ij c_i^dagger).
    Diagonalizing sum a^dagger a fixes its vacuum; occupied states follow by
    creation. This does not use any solver, normal ordering, or BdG eigensolver.
    """
    u, v = linear_coordinates(model)
    n = len(cs)
    assert_allclose(u.conj().T @ u + v.T @ v.conj(), np.eye(n), atol=2e-15, rtol=0.)
    assert_allclose(u.conj().T @ v + v.T @ u.conj(), 0., atol=2e-15, rtol=0.)
    matrices = [c.toarray() for c in cs]
    qp = [sum(u[i, j].conjugate() * matrices[i] + v[i, j] * matrices[i].conj().T
              for i in range(n)) for j in range(n)]
    count = sum(a.conj().T @ a for a in qp)
    energies, vectors = np.linalg.eigh(count)
    assert_allclose(energies[:2], [0., 1.], atol=4e-15, rtol=0.)
    unitary = np.zeros_like(count)
    unitary[:, 0] = vectors[:, 0]
    for state in range(1, 2**n):
        bit = state & -state
        mode = bit.bit_length() - 1
        unitary[:, state] = qp[mode].conj().T @ unitary[:, state ^ bit]
    assert_allclose(unitary.conj().T @ unitary, np.eye(2**n), atol=6e-15, rtol=0.)
    return unitary


def exterior_isometry(rotation):
    """Second quantization of a rectangular one-particle isometry via minors."""
    full, kept = rotation.shape
    rows, columns, values = [0], [0], [1.]
    for particles in range(1, kept + 1):
        for source in combinations(range(kept), particles):
            for target in combinations(range(full), particles):
                determinant = np.linalg.det(rotation[np.ix_(target, source)])
                if abs(determinant) > 1e-14:
                    rows.append(sum(1 << i for i in target))
                    columns.append(sum(1 << i for i in source))
                    values.append(determinant)
    return csr_matrix((values, (rows, columns)), shape=(2**full, 2**kept))


@pytest.mark.parametrize("u, detuning, field", [(1.7, .21, -.13), (-.8, -.34, .27)])
def test_anderson_atomic_spectrum_and_observables(u, detuning, field):
    impurity = Impurity.anderson(u, detuning, field)
    cs = electron_operators(2)
    diagonal = [u / 2 - detuning, field / 2, -field / 2, u / 2 + detuning]
    assert_allclose(evaluate(impurity.operator, cs).toarray(), np.diag(diagonal), atol=2e-16, rtol=2e-15)
    for name, expected in {
        "impurity_charge": [0, 1, 1, 2], "impurity_spin_z": [0, .5, -.5, 0],
        "double_occupancy_0": [0, 0, 0, 1],
    }.items():
        assert_allclose(evaluate(impurity.observables[name], cs).toarray(), np.diag(expected),
                        atol=0., rtol=0.)


def test_complex_antisymmetrized_integrals_against_pair_space_matrix():
    rng = np.random.default_rng(92026)
    n = 4
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    one_body = (raw + raw.conj().T) / 5
    pairs = list(combinations(range(n), 2))
    raw_pairs = rng.normal(size=(len(pairs), len(pairs))) + 1j * rng.normal(size=(len(pairs), len(pairs)))
    pair_hamiltonian = (raw_pairs + raw_pairs.conj().T) / 7
    interaction = np.zeros((n,) * 4, complex)
    for p, (a, b) in enumerate(pairs):
        for q, (c, d) in enumerate(pairs):
            value = pair_hamiltonian[p, q]
            interaction[a, b, c, d] = interaction[b, a, d, c] = value
            interaction[b, a, c, d] = interaction[a, b, d, c] = -value
    constant = .31
    impurity = Impurity.from_integrals(one_body, interaction, constant)
    cs = [c.toarray() for c in electron_operators(n)]
    expected = constant * np.eye(2**n, dtype=complex)
    for a, b in product(range(n), repeat=2):
        expected += one_body[a, b] * cs[a].conj().T @ cs[b]
    pair_creators = [cs[a].conj().T @ cs[b].conj().T for a, b in pairs]
    for p, q in product(range(len(pairs)), repeat=2):
        expected += pair_hamiltonian[p, q] * pair_creators[p] @ pair_creators[q].conj().T
    actual = evaluate(impurity.operator, electron_operators(n)).toarray()
    assert_allclose(actual, expected, atol=2e-15, rtol=2e-15)
    assert_allclose(actual, actual.conj().T, atol=0., rtol=0.)
    two_particle = [sum(1 << i for i in pair) for pair in pairs]
    interaction_only = Impurity.from_integrals(np.zeros((n, n)), interaction)
    matrix = evaluate(interaction_only.operator, electron_operators(n)).toarray()
    assert_allclose(matrix[np.ix_(two_particle, two_particle)], pair_hamiltonian, atol=0., rtol=0.)
    for orbital in range(2):
        diagonal = [float(bool(state & (1 << (2 * orbital))) and bool(state & (1 << (2 * orbital + 1))))
                    for state in range(16)]
        assert_allclose(evaluate(impurity.observables[f"double_occupancy_{orbital}"],
                                 electron_operators(n)).toarray(), np.diag(diagonal), atol=0., rtol=0.)


def test_one_body_integrals_have_subset_sum_spectrum_and_preserve_custom_observables():
    h = np.array([[.3, .2j], [-.2j, -.7]])
    impurity = Impurity.from_integrals(h, constant=.4)
    eigenvalues = np.linalg.eigvalsh(h)
    expected = sorted([.4, .4 + eigenvalues[0], .4 + eigenvalues[1], .4 + eigenvalues.sum()])
    assert_allclose(np.linalg.eigvalsh(evaluate(impurity.operator, electron_operators(2)).toarray()),
                    expected, atol=3e-16, rtol=2e-15)
    custom = {"impurity_charge": 7 * number(0)}
    other = Impurity(1, number(0), custom)
    custom.clear()
    assert other.observables["impurity_charge"].terms == (7 * number(0)).terms
    assert "impurity_spin_z" in other.observables


@pytest.mark.parametrize("orbitals", [0, 1.5])
def test_impurity_orbitals_validation(orbitals):
    with pytest.raises(ValueError, match="positive integer"):
        Impurity(orbitals, FermionOperator())


def test_anderson_rejects_nonfinite_parameters():
    with pytest.raises(ValueError, match="parameters must be finite"):
        Impurity.anderson(u=np.nan)


@pytest.mark.parametrize("h", [np.zeros((2, 3)), np.zeros((3, 3)), [[np.nan, 0], [0, 0]]])
def test_one_body_integral_validation(h):
    with pytest.raises(ValueError, match="finite even-dimensional square matrix"):
        Impurity.from_integrals(h)


@pytest.mark.parametrize("interaction", [np.zeros((2, 2)), np.full((2,) * 4, np.inf)])
def test_two_body_integral_validation(interaction):
    with pytest.raises(ValueError, match="integral tensor"):
        Impurity.from_integrals(np.eye(2), interaction)


def test_nonhermitian_integrals_are_rejected():
    with pytest.raises(ValueError, match="not Hermitian"):
        Impurity.from_integrals([[0., 1j], [0., 0.]])


def base_arguments():
    bath = DiscreteBath([0.], [.7], .9, 3.)
    return dict(impurity=Impurity.anderson(), baths=[bath, bath],
                tunneling=[.23 * np.eye(2), .23 * np.eye(2)])


@pytest.mark.parametrize("changes, message", [
    ({"baths": []}, "at least one DiscreteBath"),
    ({"phases": [0.]}, "finite phase"),
    ({"phase_velocities": [0., np.inf]}, "finite phase"),
    ({"tunneling": [np.zeros((2, 1)), np.eye(2)]}, "tunneling matrix"),
    ({"tunneling": [np.full((2, 2), np.nan), np.eye(2)]}, "tunneling matrix"),
    ({"direct": {(1, 0): np.eye(2)}}, "direct hopping"),
    ({"direct": {(0, 1): np.eye(3)}}, "direct hopping"),
    ({"direct_derivatives": {"phase_derivative": {}}}, "observable names"),
    ({"direct_derivatives": {"test": {(0, 1): np.eye(3)}}}, "contact observables"),
    ({"compress_pairs": True}, "compression requires"),
    ({"bath_reference": "invalid"}, "bath_reference must"),
    ({"coefficient_tolerance": -1.}, "nonnegative and finite"),
    ({"eta_basis": True, "bath_reference": "coupled"}, "isolated-reservoir QPs"),
])
def test_make_model_argument_validation(changes, message):
    arguments = base_arguments()
    arguments.update(changes)
    with pytest.raises(ValueError, match=message):
        make_model(**arguments)


@pytest.mark.parametrize("changes", [
    pytest.param({"impurity": Impurity.anderson(detuning=.1)}, id="asymmetric-gate"),
    pytest.param({"baths": [DiscreteBath([.1], [.7], .9, 3.)] * 2}, id="unpaired-band"),
    pytest.param({"tunneling": [.2 * np.eye(2), .3 * np.eye(2)]}, id="unequal-contacts"),
    pytest.param({"tunneling": [np.diag([.2, .3])] * 2}, id="spin-dependent-contact"),
    pytest.param({"direct": {(0, 1): .1 * np.eye(2)}}, id="normal-interlead-hopping"),
])
def test_eta_symmetry_requirements(changes):
    args = base_arguments()
    args.update(changes)
    with pytest.raises(ValueError, match="eta basis requires"):
        make_model(**args, eta_basis=True)


def test_compression_rejects_finite_spin_orbit_contact():
    with pytest.raises(ValueError, match="finite W_S couples"):
        make_model(**base_arguments(), eta_basis=True, compress_pairs=True,
                   direct={(0, 1): np.diag([.1j, -.1j])})


@pytest.mark.parametrize("kind", ["spin-flip", "polarized", "gapless"])
def test_coupled_vacuum_rejects_spin_mixing_polarization_and_gap_closure(kind):
    args = base_arguments()
    if kind == "spin-flip":
        args["direct"] = {(0, 1): [[0., .2j], [0., 0.]]}
        message = "spin-conserving"
    elif kind == "polarized":
        args["baths"] = [DiscreteBath([1.], [1.], .2, 3.)] * 2
        args["direct"] = {(0, 1): np.diag([2., 0.])}
        message = "unpolarized"
    else:
        args["baths"] = [DiscreteBath([0.], [1.], 1., 3.)] * 2
        args["direct"] = {(0, 1): np.diag([1. + 1e-14, -1. - 1e-14])}
        message = "nonzero excitation gap"
    with pytest.raises(ValueError, match=message):
        make_model(**args, bath_reference="coupled")


@pytest.mark.parametrize("reference", ["isolated", "coupled"])
def test_nonzero_coupled_hopping_full_matrices_spectra_and_physical_responses(reference):
    impurity = Impurity.anderson(u=1.4, detuning=.11, field=.05)
    baths = [DiscreteBath([.2], [.6], .8, 3.), DiscreteBath([-.45], [.9], .67, 2.7)]
    tunneling = [np.diag([.28 + .12j, .21 - .07j]), np.diag([.13 - .08j, .19 + .03j])]
    phases, velocities = np.array([-.6, .3]), np.array([-.4, .7])
    direct = {(0, 1): np.diag([.33 + .17j, .28 - .08j])}
    response = np.array([[.17 + .23j, .09j], [-.11, -.07 + .19j]])
    h = make_model(impurity, baths, tunneling, phases, direct, velocities,
                   bath_reference=reference, direct_derivatives={"contact_response": {(0, 1): response}})
    cs = electron_operators(6)
    matrix = electron_hamiltonian(impurity, baths, tunneling, phases, direct).toarray()
    unitary = electron_coordinate_unitary(h, cs)
    actual = evaluate(h.operator, cs).toarray()
    assert_allclose(actual, unitary.conj().T @ matrix @ unitary, atol=8e-15, rtol=3e-14)
    for parity in (0, 1):
        for sz in range(-3, 4):
            indices = sector_indices(6, parity, sz)
            if not len(indices):
                continue
            h.check_sector(Sector(parity, sz))
            assert_allclose(np.linalg.eigvalsh(actual[np.ix_(indices, indices)]),
                            np.linalg.eigvalsh(matrix[np.ix_(indices, indices)]), atol=9e-15, rtol=2e-13)
    expected_observables = {
        "impurity_charge": (cs[0].getH() @ cs[0] + cs[1].getH() @ cs[1]).toarray(),
        "impurity_spin_z": (cs[0].getH() @ cs[0] - cs[1].getH() @ cs[1]).toarray() / 2,
        "double_occupancy_0": (cs[0].getH() @ cs[0] @ cs[1].getH() @ cs[1]).toarray(),
    }
    step = 2e-6
    expected_observables["phase_derivative"] = (
        electron_hamiltonian(impurity, baths, tunneling, phases + step * velocities, direct)
        - electron_hamiltonian(impurity, baths, tunneling, phases - step * velocities, direct)
    ).toarray() / (2 * step)
    contact = np.zeros_like(matrix)
    for s, t in product(range(2), repeat=2):
        term = (np.exp(.5j * (phases[0] - phases[1])) * response[s, t]
                * np.sqrt(baths[0].weights[0] * baths[1].weights[0])
                * (cs[2 + s].getH() @ cs[4 + t]).toarray())
        contact += term + term.conj().T
    expected_observables["contact_response"] = contact
    for name, observable in expected_observables.items():
        assert_allclose(evaluate(h.observables[name], cs).toarray(),
                        unitary.conj().T @ observable @ unitary,
                        atol=6e-11 if name == "phase_derivative" else 5e-15,
                        rtol=2e-10 if name == "phase_derivative" else 0.)
    # An arbitrary electron polynomial probes the inverse map beyond bilinears.
    physical = .2j * create(0) * annihilate(3) * create(5) + .13 * annihilate(4) + .1j
    assert_allclose(evaluate(h.physical_operator(physical), cs).toarray(),
                    unitary.conj().T @ evaluate(physical, cs).toarray() @ unitary, atol=1e-15, rtol=2e-14)
    if reference == "coupled":
        # Coupling genuinely changes the reservoir vacuum; this is not the W=0 path.
        assert abs(h.metadata["bath_vacuum_shift"]) > 1e-3
        bath_impurity = Impurity(1, FermionOperator())
        decoupled = make_model(bath_impurity, baths, [np.zeros((2, 2))] * 2,
                               phases, direct, bath_reference="coupled")
        assert all(not word or (len(word) == 2 and word[0] == -word[1])
                   for word in decoupled.operator.terms)
        bath_matrix = electron_hamiltonian(bath_impurity, baths, [np.zeros((2, 2))] * 2,
                                           phases, direct).toarray()
        assert_allclose(decoupled.operator.terms[()], np.linalg.eigvalsh(bath_matrix)[0],
                        atol=4e-15, rtol=2e-13)


@pytest.mark.slow
@pytest.mark.convergence
def test_unequal_coupled_baths_against_electron_sectors():
    """A fixed eight-mode problem exercises rectangular interlead contact blocks."""
    impurity = Impurity.anderson(u=1.6, detuning=.11, field=.09)
    baths = [DiscreteBath([-.8, .6], [.3, .4], .95, 3.), DiscreteBath([.2], [.75], .8, 2.5)]
    tunneling = [np.diag([.26 + .08j, .21 - .05j]), np.diag([.18 - .06j, .3 + .03j])]
    phases = [-.6, .9]
    direct = {(0, 1): np.diag([.24 + .12j, .19 - .08j])}
    cs = electron_operators(8)
    physical = electron_hamiltonian(impurity, baths, tunneling, phases, direct).toarray()
    physical_contact = physical - electron_hamiltonian(impurity, baths, tunneling, phases).toarray()
    assert np.linalg.norm(physical_contact) > .1
    physical_charge = (cs[0].getH() @ cs[0] + cs[1].getH() @ cs[1]).toarray()

    def thermal_sector(matrix, observables, indices):
        # Thermal traces avoid matching arbitrary eigenvectors in degeneracies.
        energies, vectors = np.linalg.eigh(matrix[np.ix_(indices, indices)])
        probabilities = np.exp(-1.3 * (energies - energies[0]))
        probabilities /= probabilities.sum()
        means = [probabilities @ np.diag(vectors.conj().T @ obs[np.ix_(indices, indices)] @ vectors)
                 for obs in observables]
        return energies, means

    for reference in ("isolated", "coupled"):
        h = make_model(impurity, baths, tunneling, phases, direct, bath_reference=reference,
                       direct_derivatives={"contact_hopping": direct})
        actual = evaluate(h.operator, cs).toarray()
        observables = [evaluate(h.observables[name], cs).toarray()
                       for name in ("impurity_charge", "contact_hopping")]
        for sz in range(-4, 5):
            indices = sector_indices(8, sz % 2, sz)
            context = f"reference={reference}, twice_sz={sz}"
            h.check_sector(Sector(sz % 2, sz))
            expected = thermal_sector(physical, [physical_charge, physical_contact], indices)
            result = thermal_sector(actual, observables, indices)
            for computed, independent in zip(result, expected, strict=True):
                assert_allclose(computed, independent, atol=2e-13, rtol=2e-13, err_msg=context)


def test_eta_coordinates_with_finite_spin_orbit_contact_against_electron_matrix():
    bath = DiscreteBath([0.], [.7], .9, 3.)
    impurity = Impurity.anderson(u=1.7, field=.13)
    hopping, phases = [.24 * np.eye(2)] * 2, [.31, 1.04]
    direct = {(0, 1): np.diag([.19j, -.19j])}
    h = make_model(impurity, [bath] * 2, hopping, phases, direct, eta_basis=True)
    cs = electron_operators(6)
    unitary = electron_coordinate_unitary(h, cs)
    physical = electron_hamiltonian(impurity, [bath] * 2, hopping, phases, direct).toarray()
    actual = evaluate(h.operator, cs).toarray()
    assert_allclose(actual, unitary.conj().T @ physical @ unitary, atol=2e-14, rtol=2e-14)
    charge = (cs[0].getH() @ cs[0] + cs[1].getH() @ cs[1]).toarray()
    assert_allclose(evaluate(h.observables["impurity_charge"], cs).toarray(),
                    unitary.conj().T @ charge @ unitary, atol=2e-14, rtol=2e-14)
    occupations = np.array([[(state >> i) & 1 for i in range(6)] for state in range(64)])
    charges = [occupations.sum(axis=1) % 2, occupations @ h.spins,
               np.prod(np.where(occupations, h.eta_labels, 1), axis=1)]
    for charge in charges:
        assert_allclose(actual[charge[:, None] != charge[None, :]], 0., atol=1e-14, rtol=0.)
    h.check_sector(Sector(1, 1, 1))


@pytest.mark.parametrize("xi, phi", [([0.], 0.), ([-.5, .5], .83)])
def test_compressed_dark_vacuum_projection_including_contact_contractions(xi, phi):
    bath = DiscreteBath(xi, [.6] * len(xi), .9, 3.)
    impurity = Impurity.anderson(u=1.3, field=.07)
    phases = [.27 - phi / 2, .27 + phi / 2]
    responses = {
        "normal_contact": {(0, 1): np.eye(2)},
        "spin_orbit_contact": {(0, 1): np.diag([1j, -1j])},
        "spin_mixed_contact": {(0, 1): np.array([[.2 + .3j, .4j], [-.1, -.3 + .2j]])},
    }
    arguments = dict(impurity=impurity, baths=[bath] * 2, tunneling=[.21 * np.eye(2)] * 2,
                     phases=phases, phase_velocities=[-.5, .5], eta_basis=True,
                     direct_derivatives=responses)
    full = make_model(**arguments)
    reduced = make_model(**arguments, compress_pairs=True)
    uf, vf = linear_coordinates(full)
    ur, vr = linear_coordinates(reduced)
    rotation = uf.conj().T @ ur + vf.T @ vr.conj()
    assert_allclose(uf.conj().T @ vr + vf.T @ ur.conj(), 0., atol=5e-16, rtol=0.)
    assert_allclose(rotation.conj().T @ rotation, np.eye(len(reduced.spins)), atol=1e-15, rtol=0.)
    isometry = exterior_isometry(rotation)
    assert_allclose((isometry.getH() @ isometry).toarray(), np.eye(isometry.shape[1]), atol=2e-14, rtol=0.)
    full_cs, reduced_cs = electron_operators(len(full.spins)), electron_operators(len(reduced.spins))
    operators = [(full.operator, reduced.operator)] + [
        (full.observables[name], value) for name, value in reduced.observables.items()]
    for original, compressed in operators:
        projected = isometry.getH() @ evaluate(original, full_cs) @ isometry
        assert_allclose(evaluate(compressed, reduced_cs).toarray(), projected.toarray(),
                        atol=9e-15, rtol=2e-14)
    dark = full.bath_modes - reduced.bath_modes
    assert dark > 0 and len(reduced.metadata["decoupled_qp_energies"]) == dark
    assert_allclose(reduced.metadata["decoupled_qp_energies"], bath.energies[0], atol=0., rtol=0.)
    # Products of separately projected electrons really do miss a finite
    # constant; equality above therefore exercises the dark-vacuum correction.
    maps = [FermionOperator.from_records(record) for record in
            reduced.metadata["coordinates"]["electron_annihilators"]]
    left = [sum(np.sqrt(w) * maps[2 + 2 * i + s] for i, w in enumerate(bath.weights)) for s in range(2)]
    right = [sum(np.sqrt(w) * maps[2 + 2 * bath.levels + 2 * i + s]
                 for i, w in enumerate(bath.weights)) for s in range(2)]
    naive = sum((np.exp(-.5j * phi) * left[s].dagger() * right[s] for s in range(2)), FermionOperator())
    naive = naive + naive.dagger()
    missed = evaluate(reduced.observables["normal_contact"] - naive, reduced_cs).toarray()
    assert abs(np.trace(missed) / len(missed)) > .01
    assert_allclose(missed, np.eye(len(missed)) * missed[0, 0], atol=2e-15, rtol=0.)
    with pytest.raises(ValueError, match="uncompressed"):
        reduced.physical_operator(number(0))


@pytest.mark.parametrize("kwargs, message", [
    ({"geometry": "chain"}, "geometry must"), ({"gamma": -.1}, "Gamma nonnegative"),
    ({"phi": np.nan}, "parameters must be finite"),
    ({"geometry": "single", "rho_ws": .1}, "two reservoirs"),
])
def test_reference_model_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        reference_model(DiscreteBath([0.], [1.]), **kwargs)


@pytest.mark.parametrize("geometry, symmetry, rho_ws, rho_wn, detuning", [
    ("single", True, 0., 0., .13), ("junction", False, .02, .03, .11),
    ("junction", True, .03, 0., 0.), ("junction", True, 0., 0., 0.),
])
def test_reference_model_normalizations_and_symmetry_choices(geometry, symmetry, rho_ws, rho_wn, detuning):
    bath = DiscreteBath([0.], [.7], .9, 3.)
    gamma, phi = .31, .81
    h = reference_model(bath, gamma=gamma, phi=phi, rho_ws=rho_ws, rho_wn=rho_wn,
                        geometry=geometry, symmetry=symmetry, detuning=detuning)
    leads = 1 if geometry == "single" else 2
    encoded = np.asarray(h.metadata["tunneling"])
    hopping = encoded[..., 0] + 1j * encoded[..., 1]
    assert_allclose(sum(np.pi * bath.rho * np.diag(t.conj().T @ t) for t in hopping),
                    gamma * np.ones(2), atol=5e-17, rtol=2e-16)
    assert h.metadata["geometry"] == ("single-reservoir" if leads == 1 else "2-reservoir")
    eta_expected = leads == 2 and symmetry and not rho_wn and not detuning
    assert (h.eta_labels is not None) == eta_expected
    assert h.metadata["paired_mode_compression"] == bool(eta_expected and not rho_ws)
    assert ("phase_derivative" in h.observables) == (leads == 2)
    assert ("rho_ws_derivative" in h.observables) == (leads == 2)
    assert h.metadata["phi"] == (phi if leads == 2 else 0.)
