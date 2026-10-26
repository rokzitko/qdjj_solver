"""Probe matrix elements between static and Floquet states, spectral lines and line shapes."""

import dataclasses

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.special import jv

from qdjj_solver import (FermionOperator, FloquetTransitions, Hamiltonian, Impurity, ParameterDrive,
                         PeriodicHamiltonian, Sector, SpectralLines, annihilate, cosh_grid, create,
                         hermitian_pair, line_shape, number, reference_model)
from qdjj_solver.common.problem import check_transition
from qdjj_solver.common.spectra import degenerate_groups
from qdjj_solver.qp_solver import ProjectedOperator, SolverOptions, floquet_propagator, solve, solve_floquet
from electron_oracle import electron_hamiltonian, electron_operators, sector_indices


def junction(**changes):
    base = dict(bath=cosh_grid(pairs=1, bandwidth=3.), u=1., gamma=.5, phi=.8, detuning=.2, field=0.,
                symmetry=False, compress=False)
    base.update(changes)
    return base


def drive(amplitude_phi=.6, amplitude_eps=.2, theta=.3, omega=.9, **changes):
    drives = (ParameterDrive.cosine("phi", amplitude_phi), ParameterDrive.cosine("detuning", amplitude_eps, theta))
    return PeriodicHamiltonian.from_parameters(reference_model, junction(**changes), drives, omega, samples=32)


def lines_as_set(lines, digits=9):
    return sorted((round(f, digits), round(s, digits)) for f, s in zip(lines.frequencies, lines.strengths,
                                                                       strict=True))


def test_static_matrix_elements_match_electron_oracle():
    model = junction()
    h = reference_model(**model)
    even = solve(h, None, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2))
    odd = solve(h, None, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1))
    qp = abs(odd.matrix_elements("d_up_dagger", even))**2
    bath = model["bath"]
    t = np.sqrt(model["gamma"]/(np.pi*bath.rho*2))*np.eye(2)
    matrix = electron_hamiltonian(Impurity.anderson(1., .2), [bath, bath], [t, t], [-.4, .4])
    cs = electron_operators(10)
    ie, io = sector_indices(10, 0, 0), sector_indices(10, 1, 1)
    _, ve = np.linalg.eigh(matrix[ie][:, ie].toarray())
    _, vo = np.linalg.eigh(matrix[io][:, io].toarray())
    oracle = abs(vo[:, :1].conj().T @ cs[0].getH()[io][:, ie].toarray() @ ve[:, :2])**2
    assert_allclose(qp, oracle, rtol=0, atol=1e-12)
    # Same-sector elements agree with ProjectedOperator, and overlaps are orthonormal.
    current = ProjectedOperator(h.observables["phase_derivative"], even.basis)
    assert_allclose(even.matrix_elements("current")[0, 1], current.expectation(even.vectors[:, 0], even.vectors[:, 1]),
                    atol=1e-14)
    assert_allclose(even.overlaps(), np.eye(2), atol=1e-12)
    truncated = solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2))
    assert abs(truncated.overlaps(even)[0, 0]) > .99
    with pytest.raises(ValueError, match="no term"):
        odd.matrix_elements("d_up", even)
    other = solve(reference_model(**{**model, "symmetry": True, "detuning": 0.}), 2, Sector(0, 0))
    with pytest.raises(ValueError, match="mode labels"):
        even.matrix_elements("charge", other)
    with pytest.raises(ValueError, match="QP Eigenstates"):
        even.matrix_elements("charge", object())


def test_probe_resolution():
    h = reference_model(**junction())
    assert h.probe("current") is h.observables["phase_derivative"]
    charge = h.probe("charge", 0) - h.probe("charge")
    assert max(map(abs, charge.terms.values()), default=0.) < 1e-14
    spin_x = h.probe("spin_x")
    assert spin_x.hermiticity_error() < 1e-14
    assert max(map(abs, (h.probe("spin_plus") + h.probe("spin_minus") - spin_x*2).terms.values()), default=0.) < 1e-14
    assert h.probe("double_occupancy_0") is h.observables["double_occupancy_0"]
    for name, orbital, message in (("unknown", None, "unknown probe"), ("current", 0, "not of an orbital"),
                                   ("charge", 1, "orbital"), (3, None, "string")):
        with pytest.raises(ValueError, match=message):
            h.probe(name, orbital)
    single = reference_model(cosh_grid(pairs=1, bandwidth=3.), geometry="single", compress=False, symmetry=False)
    with pytest.raises(ValueError, match="phase_velocities"):
        single.probe("current")
    with pytest.raises(ValueError, match="uncompressed"):
        reference_model(cosh_grid(pairs=1, bandwidth=3.), detuning=0.).probe("d_up")
    periodic = drive()
    assert set(periodic.probe("current")) >= {-1, 0, 1}
    assert set(periodic.probe("charge")) == {0}
    assert set(periodic.probe("d_up")) == {0}
    bare = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", .5), 1.,
                                               observables=False)
    with pytest.raises(ValueError, match="not sampled"):
        bare.probe("current")
    with pytest.raises(ValueError, match="no term"):
        check_transition(h.probe("d_up"), h.spins, None, Sector(0, 0), Sector(1, 1))
    with pytest.raises(ValueError, match="zero"):
        check_transition(FermionOperator(), h.spins, None, Sector(0, 0), Sector(0, 0))


def _propagate_exact(model, amplitudes, theta, omega, sector, times, start):
    """Independent propagation with H(t) rebuilt from the model at every Gauss point."""
    def hamiltonian(t):
        h = reference_model(**{**model, "phi": model["phi"] + amplitudes[0]*np.cos(omega*t),
                               "detuning": model["detuning"] + amplitudes[1]*np.cos(omega*t + theta)})
        return h
    basis = start[1]
    def matrix(t):
        return ProjectedOperator(hamiltonian(t).operator, basis).sparse_matrix().toarray()
    def exponential(m, h):
        values, vectors = np.linalg.eigh((m + m.conj().T)/2)
        return (vectors*np.exp(-1j*h*values)) @ vectors.conj().T
    state, t, steps = start[0].copy(), 0., 200
    results = []
    nodes = .5 + np.array([-1., 1.])*np.sqrt(3)/6
    for target in times:
        h = (target - t)/steps
        for _ in range(steps):
            h1, h2 = (matrix(t + c*h) for c in nodes)
            generator = (h1 + h2)/2 + 1j*np.sqrt(3)/12*h*(h1 @ h2 - h2 @ h1)
            state = exponential(generator, h) @ state
            t += h
        current = ProjectedOperator(hamiltonian(t).observables["phase_derivative"], basis).sparse_matrix()
        results.append(state.conj().T @ (current @ state))
    return results


def test_elements_match_time_propagation():
    amplitudes, theta, omega = (.6, .2), .3, .9
    periodic = drive(*amplitudes, theta, omega)
    sector = Sector(0, 0)
    states = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=10, selection="all", method="dense",
                                                     dense_dimension=4000))
    pick = np.array([states.dressed_state(0), states.dressed_state(1)])
    transitions = states.matrix_elements("current", bra_states=pick, ket_states=pick)
    times = [.37*periodic.period, 1.21*periodic.period]
    exact = _propagate_exact(junction(), amplitudes, theta, omega, sector, times,
                             (states.state_at(0.)[:, pick], states.basis))
    eps = states.quasienergies[pick]
    for t, direct in zip(times, exact, strict=True):
        series = np.einsum("q,qab->ab", np.exp(1j*transitions.photons*omega*t), transitions.elements)
        series *= np.exp(1j*(eps[:, None] - eps[None, :])*t)
        assert_allclose(series, direct, rtol=0, atol=1e-8)
        wrong = np.einsum("q,qab->ab", np.exp(-1j*transitions.photons*omega*t), transitions.elements)
        wrong *= np.exp(1j*(eps[:, None] - eps[None, :])*t)
        assert np.max(abs(wrong - direct)) > 1e-3


def test_sum_rule_hermiticity_and_replicas():
    periodic = drive()
    sector = Sector(0, 0)
    states = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=10, selection="all", method="dense",
                                                     dense_dimension=4000))
    assert states.metadata["complete"]
    current = states.matrix_elements("current")
    assert current.complete
    assert_allclose(current.completeness, 1., atol=1e-10)
    # Hermitian, time-dependent probe: conj(X^(q)_ab) = X^(-q)_ba.
    flipped = current.elements[::-1].transpose(0, 2, 1)
    assert_allclose(np.conj(current.elements), flipped, atol=1e-12)
    assert_allclose(current.photons, -current.photons[::-1])
    # The bra replica (eps+omega, u_{m-1}) has X'^(q) = X^(q+1) and the same lines.
    shifted = dataclasses.replace(states, quasienergies=states.quasienergies + periodic.omega,
                                  harmonic_indices=states.harmonic_indices + 1)
    moved = shifted.matrix_elements("current", states)
    assert_allclose(moved.element(2), current.element(3), atol=1e-12)
    def table(transitions):
        keep = transitions.strengths > 1e-14
        return sorted(zip(np.round(transitions.frequencies[keep], 9), np.round(transitions.strengths[keep], 9), strict=True))
    assert table(moved) == pytest.approx(table(current), abs=1e-9)
    # Tracked subsets capture less than the complete set.
    tracked = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=10, states=2))
    partial = tracked.matrix_elements("current")
    assert not partial.complete and np.all(partial.completeness < 1)


def test_cross_sector_tien_gordon_both_gauges():
    amplitude, omega, theta = .7, 1., .3
    model = junction(gamma=0.)
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                   ParameterDrive.cosine("detuning", amplitude, theta), omega,
                                                   samples=16)
    for drop in (False, True):
        options = dict(floquet_cutoff=12, selection="all", method="dense", drop_scalar_harmonics=drop,
                       dense_dimension=4000)
        even = solve_floquet(periodic, 2, Sector(0, 0), options)
        odd = solve_floquet(periodic, 2, Sector(1, 1), options)
        addition = odd.matrix_elements("d_up_dagger", even)
        removal = even.matrix_elements("d_up", odd)
        assert_allclose(addition.completeness, 1., atol=1e-10)
        # Both static states with charge 0 and 1 couple; each line carries J_q(A/omega)^2.
        lines = addition.lines(initial=even.dressed_state(0))
        weights = {int(q): s for q, s in zip(lines.photons, lines.strengths, strict=True)}
        total = sum(weights.values())
        for q in (-2, -1, 0, 1, 2):
            assert_allclose(weights.get(q, 0.)/total, jv(q, amplitude/omega)**2, atol=1e-10)
        # d and d^dagger: X^(q)_{odd,even}(d^dagger) = conj((d)^(-q)_{even,odd}).
        assert_allclose(addition.elements, np.conj(removal.elements[::-1].transpose(0, 2, 1)), atol=1e-12)
        coupled = abs(addition.element(0)) > 1e-6
        ratio, reference = addition.element(1)[coupled], addition.element(0)[coupled]
        assert_allclose(ratio/reference, jv(1, amplitude/omega)/jv(0, amplitude/omega)*np.exp(1j*theta), atol=1e-10)
    # Bra and ket computed in different scalar gauges would give silently wrong elements.
    kept = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=12, selection="all", method="dense"))
    with pytest.raises(ValueError, match="scalar gauges"):
        odd.matrix_elements("d_up_dagger", kept)


def test_adiabatic_gate_drive_follows_generalized_tien_gordon():
    # Omega = 0.05 lies far below the gaps within each parity sector (0.66 and 1.51), so the
    # D <- S tunnelling sidebands follow the Fourier series of M(eps(t)) exp[i int dE(eps(t'))dt']
    # built from static calculations; the residual difference is an O(Omega/gap) correction.
    amplitude, omega = .2, .05
    model = junction()
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                   ParameterDrive.cosine("detuning", amplitude), omega, samples=16)
    harmonics = int(np.ceil(amplitude/omega)) + 10
    even = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=harmonics, states=1))
    odd = solve_floquet(periodic, 2, Sector(1, 1), dict(floquet_cutoff=harmonics, selection="all", method="dense",
                                                        dense_dimension=8000))
    table = odd.matrix_elements("d_up_dagger", even)
    assert_allclose(table.completeness, 1., atol=1e-8)
    lines = table.lines(initial=even.dressed_state(0))
    doublet = odd.dressed_state(0)
    floquet = {}
    for q, strength, final in zip(lines.photons, lines.strengths, lines.final, strict=True):
        if doublet in final:
            floquet[int(q)] = floquet.get(int(q), 0.) + strength
    samples = 128
    theta = 2*np.pi*np.arange(samples)/samples
    gap, element = [], []
    for value in model["detuning"] + amplitude*np.cos(theta):
        h = reference_model(**{**model, "detuning": float(value)})
        s = solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=1))
        d = solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1))
        gap.append(d.energies[0] - s.energies[0])
        element.append(abs(d.matrix_elements("d_up_dagger", s)[0, 0]))
    # Spectral integration of the dynamical phase, dt = dtheta/omega.
    n = np.fft.fftfreq(samples, 1/samples)
    coefficients = np.fft.fft(np.asarray(gap) - np.mean(gap))
    phase = np.real(np.fft.ifft(np.divide(coefficients, 1j*n, out=np.zeros(samples, complex), where=n != 0)))/omega
    adiabatic = np.fft.fft(np.asarray(element)*np.exp(1j*phase))/samples
    static = solve(reference_model(**model), 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=1))
    static_odd = solve(reference_model(**model), 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1))
    norm = abs(static_odd.matrix_elements("d_up_dagger", static)[0, 0])**2
    for q in (-2, -1, 0, 1, 2):
        assert_allclose(floquet.get(q, 0.)/norm, abs(adiabatic[q % samples])**2/norm, atol=.012)
    # Plain Tien-Gordon is symmetric in q; the junction sidebands are not.
    assert floquet[1] - floquet[-1] > .04*norm


def test_circular_two_level_lines():
    for d0, amplitude, omega in ((1., .3, 1.2), (1., .3, .8), (1., .05, 1.)):
        def model(t, d0=d0, amplitude=amplitude, omega=omega):
            drive_term = amplitude/2*complex(np.exp(-1j*omega*t))*create(0)*annihilate(1)
            return Hamiltonian(d0/2*(number(0) - number(1)) + hermitian_pair(drive_term), 2, (1, -1))
        periodic = PeriodicHamiltonian.from_callable(model, omega, samples=8, assume_fixed_coordinates=True)
        states = solve_floquet(periodic, 0, Sector(1, None), dict(floquet_cutoff=6, selection="all", method="dense"))
        sigma_x = hermitian_pair(create(0)*annihilate(1))
        table = states.matrix_elements(sigma_x)
        width = np.hypot(d0 - omega, amplitude)
        angle = np.arctan2(amplitude, d0 - omega)
        amplitudes = {+1: (np.cos(angle/2), np.sin(angle/2)), -1: (-np.sin(angle/2), np.cos(angle/2))}
        branch = {+1: (omega + width)/2, -1: (omega - width)/2}
        expected = []
        for a in (+1, -1):
            for b in (+1, -1):
                up_a, down_a = amplitudes[a]
                up_b, down_b = amplitudes[b]
                for q, value in ((1, up_a*down_b), (-1, down_a*up_b)):
                    expected.append((branch[a] - branch[b] + q*omega, abs(value)**2))
        computed = []
        for i in range(len(table.photons)):
            for a in range(2):
                for b in range(2):
                    if table.strengths[i, a, b] > 1e-20:
                        computed.append((table.frequencies[i, a, b], table.strengths[i, a, b]))
        # The branch of each state is fixed by the solver; compare modulo omega and by strength.
        expected = sorted((round(np.mod(f, omega), 9), round(s, 9)) for f, s in expected)
        computed = sorted((round(np.mod(f, omega), 9), round(s, 9)) for f, s in computed)
        assert computed == pytest.approx(expected, abs=1e-9)


def test_undriven_limit_and_degenerate_merging():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", 0.), .9,
                                                   samples=8)
    sector = Sector(1, None)
    states = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=2, selection="all", method="dense"))
    table = states.matrix_elements("spin_x")
    assert np.max(abs(table.elements[table.photons != 0])) < 1e-12
    static = solve(periodic.reference, 2, sector, SolverOptions(method="dense", eigenpairs=states.basis.dimension))
    lines = table.lines(initial=states.dressed_state(0), elastic=False)
    reference = static.matrix_elements("spin_x")
    groups = degenerate_groups(static.energies)
    start = np.flatnonzero(groups == groups[0])
    expected = {}
    for a in range(len(static.energies)):
        key = groups[a]
        expected[key] = expected.get(key, 0.) + np.sum(abs(reference[a, start])**2)/len(start)
    expected.pop(groups[0], None)
    assert_allclose(sorted(lines.strengths), sorted(v for v in expected.values() if v > 1e-14), atol=1e-10)
    # Random unitary rotations inside degenerate multiplets leave the merged lines unchanged.
    rng = np.random.default_rng(5)
    modes = states.modes.copy()
    for group in np.unique(degenerate_groups(states.quasienergies, states.omega)):
        members = np.flatnonzero(degenerate_groups(states.quasienergies, states.omega) == group)
        if len(members) > 1:
            q, _ = np.linalg.qr(rng.normal(size=(len(members),)*2) + 1j*rng.normal(size=(len(members),)*2))
            modes[:, :, members] = np.einsum("mdn,nk->mdk", modes[:, :, members], q)
    rotated = dataclasses.replace(states, modes=modes).matrix_elements("spin_x")
    first = table.lines(initial=states.dressed_state(0))
    second = rotated.lines(initial=states.dressed_state(0))
    assert lines_as_set(first) == pytest.approx(lines_as_set(second), abs=1e-9)


def test_sambe_propagator_and_static_paths_agree():
    periodic = drive()
    sector = Sector(0, 0)
    options = dict(floquet_cutoff=8, selection="all", dense_dimension=4000)
    sambe = solve_floquet(periodic, 2, sector, {**options, "method": "dense"})
    stepped = floquet_propagator(periodic, 2, sector, {**options, "steps": 128})
    static = solve_floquet(periodic, 2, sector, {**options, "method": "dense", "basis": "eigen", "basis_states": 42})
    grid = np.linspace(-2.5, 2.5, 401)
    spectra = [result.matrix_elements("charge").spectrum(grid, .02, initial=result.dressed_state(0))
               for result in (sambe, stepped, static)]
    assert_allclose(spectra[1], spectra[0], rtol=0, atol=1e-7*np.max(spectra[0]))
    assert_allclose(spectra[2], spectra[0], rtol=0, atol=1e-9*np.max(spectra[0]))


def test_line_shapes_kinds_and_guards():
    x = np.linspace(-200, 200, 400001)
    for shape in ("lorentzian", "gaussian"):
        assert line_shape(0.5, .5, shape)/line_shape(0., .5, shape) == pytest.approx(.5)
    assert np.trapezoid(line_shape(x, .5, "gaussian"), x) == pytest.approx(1., abs=1e-9)
    assert np.trapezoid(line_shape(x, .5, "lorentzian"), x) == pytest.approx(1., abs=2e-3)
    with pytest.raises(ValueError, match="width"):
        line_shape(x, 0.)
    with pytest.raises(ValueError, match="shape"):
        line_shape(x, 1., "voigt")
    lines = SpectralLines(np.array([-1., .5, 1.]), np.array([.2, .3, .5]), [(0,), (1,), (2,)], [(0,)]*3,
                          np.array([0, 0, 1]), np.array([False, False, True]))
    grid = np.linspace(-3, 3, 7)
    signed = lines.spectrum(grid, .1)
    assert_allclose(lines.spectrum(grid, .1, kind="rate"), signed + signed[::-1])
    assert_allclose(lines.spectrum(grid, .1, kind="net"), signed - signed[::-1])
    with pytest.raises(ValueError, match="kind"):
        lines.spectrum(grid, .1, kind="absorption")
    assert lines.record()["photons"] == [0, 0, 1]
    assert_allclose(degenerate_groups([0., 1e-10, .5, 1. - 1e-10], omega=1.), [0, 0, 1, 0])
    periodic = drive()
    states = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=4, states=2))
    other = solve_floquet(drive(omega=1.), 2, Sector(0, 0), dict(floquet_cutoff=4, states=2))
    with pytest.raises(ValueError, match="different drives"):
        states.matrix_elements("charge", other)
    with pytest.raises(ValueError, match="no harmonic"):
        states.matrix_elements("d_up")
    with pytest.raises(ValueError, match="probe"):
        states.matrix_elements(3.)
    with pytest.raises(ValueError, match="orbital"):
        states.matrix_elements(FermionOperator.scalar(1.), orbital=0)
    with pytest.raises(ValueError, match="zero"):
        states.matrix_elements({1: FermionOperator()})
    with pytest.raises(ValueError, match="distinct"):
        states.matrix_elements("charge", ket_states=[0, 0])
    with pytest.raises(MemoryError):
        states.matrix_elements("charge", max_memory_gib=1e-9)
    table = states.matrix_elements("charge")
    assert isinstance(table, FloquetTransitions)
    with pytest.raises(ValueError, match="initial state"):
        table.lines(initial=7)
    with pytest.raises(ValueError, match="weight"):
        table.lines(initial=[0, 1], weights=[1.])
    with pytest.raises(ValueError, match="dense"):
        solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=40, selection="all", method="shift-invert"))


def test_sampled_contraction_matches_direct_sum():
    """The exact-grid DFT contraction equals the explicit sum over harmonics and probe orders."""
    periodic = drive()
    even = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=5, states=3))
    odd = solve_floquet(periodic, 2, Sector(1, 1), dict(floquet_cutoff=5, states=2, basis="eigen", basis_states=9))
    for bra, ket, probe in ((even, even, "current"), (odd, even, "d_up_dagger")):
        table = bra.matrix_elements(probe, ket)
        harmonics = periodic.probe(probe)
        u_bra, u_ket = bra.fock_modes(), ket.fock_modes()
        direct = {}
        for k, operator in harmonics.items():
            projected = ProjectedOperator(operator, ket.basis)
            for j, m_ket in enumerate(ket.harmonic_indices):
                image = np.column_stack([projected.action_between(np.ascontiguousarray(v), bra.basis)
                                         for v in u_ket[j].T])
                for i, m_bra in enumerate(bra.harmonic_indices):
                    q = m_ket + k - m_bra
                    direct[q] = direct.get(q, 0.) + u_bra[i].conj().T @ image
        for q, value in direct.items():
            assert_allclose(table.element(q), value, rtol=0, atol=1e-13)


def test_transition_checks_and_guards_for_symmetry_labels():
    h = reference_model(**junction(detuning=0., symmetry=True))
    assert h.eta_labels is not None
    hopping = h.probe("current")
    check_transition(hopping, h.spins, h.eta_labels, Sector(0, 0, 1), Sector(0, 0, 1))
    with pytest.raises(ValueError, match="no term"):
        check_transition(hopping, h.spins, h.eta_labels, Sector(0, 0, 1), Sector(0, 0, -1))
    with pytest.raises(ValueError, match="eta reduction"):
        check_transition(hopping, h.spins, None, Sector(0, 0, 1), Sector(0, 0, 1))
    pair = create(0)*create(1)
    check_transition(pair, (1, -1), None, Sector(0, 0, particle_number=0), Sector(0, 0, particle_number=2))
    with pytest.raises(ValueError, match="no term"):
        check_transition(pair, (1, -1), None, Sector(0, 0, particle_number=0), Sector(0, 0, particle_number=0))
    with pytest.raises(ValueError, match="one-dimensional"):
        degenerate_groups(np.zeros((2, 2)))
    assert len(degenerate_groups([])) == 0
    periodic = drive()
    table = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=3, states=2)).matrix_elements("charge")
    split = dataclasses.replace(table, metadata={**table.metadata, "multiplet_split": True})
    with pytest.raises(ValueError, match="multiplet"):
        split.lines()
    assert len(table.lines(merge_degenerate=False, elastic=False)) >= 1
