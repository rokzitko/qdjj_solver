"""Floquet quasienergies: drive sampling, Sambe, propagator, van Vleck and RWA checks."""

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.optimize import linear_sum_assignment
from scipy.special import jv

from qdjj_solver import (FermionOperator, Hamiltonian, Harmonic, Impurity, ParameterDrive,
                         PeriodicHamiltonian, Sector, annihilate, cosh_grid, create, hermitian_pair,
                         make_model, number, reference_model, save_result)
from qdjj_solver.common.floquet import fold_quasienergies, quasienergy_distance
from qdjj_solver.qp_solver import (FloquetOptions, ProjectedOperator, SolverOptions, floquet_propagator,
                                   rwa_model, solve, solve_floquet)
from electron_oracle import electron_hamiltonian, electron_operators, evaluate, sector_indices


def junction(**changes):
    base = dict(bath=cosh_grid(pairs=1, bandwidth=3.), u=1.6, gamma=.4, phi=.8, detuning=.1, field=0.,
                symmetry=False, compress=False)
    base.update(changes)
    return base


def single(phase=0., detuning=.15, u=1.3, gamma=.35, bandwidth=2.5):
    bath = cosh_grid(pairs=1, bandwidth=bandwidth)
    t = np.sqrt(gamma/(np.pi*bath.rho))*np.eye(2)
    return make_model(Impurity.anderson(u, detuning), [bath], [t], phases=[phase])


def difference(a, b):
    return max(map(abs, (a-b).terms.values()), default=0.)


def match(a, b, omega):
    """Largest circular distance after optimal assignment of two quasienergy sets."""
    distance = quasienergy_distance(np.asarray(a)[:, None], np.asarray(b)[None, :], omega)
    rows, columns = linear_sum_assignment(distance)
    return distance[rows, columns].max()


def circular_two_level(d0=1., amplitude=.3, omega=1.2):
    def model(t):
        drive = amplitude/2*complex(np.exp(-1j*omega*t))*create(0)*annihilate(1)
        return Hamiltonian(d0/2*(number(0)-number(1)) + hermitian_pair(drive), 2, (1, -1))
    return PeriodicHamiltonian.from_callable(model, omega, samples=8, assume_fixed_coordinates=True)


def test_parameter_drive_validation_and_values():
    drive = ParameterDrive("phi", (Harmonic(1, .3, .2), Harmonic(3, .1)))
    t, omega = .37, 1.3
    assert_allclose(drive.value(t, omega), .3*np.cos(omega*t+.2) + .1*np.cos(3*omega*t), rtol=0, atol=1e-15)
    assert drive.max_order == 3
    assert drive.record()["harmonics"][1] == dict(order=3, amplitude=.1, phase=0.)
    assert ParameterDrive.cosine("u", .2, order=2).harmonics == (Harmonic(2, .2),)
    for arguments in ((0, 1.), (1.5, 1.), (True, 1.), (1, np.inf), (1, 1j), (1, 1., np.nan)):
        with pytest.raises(ValueError, match="harmonic"):
            Harmonic(*arguments)
    with pytest.raises(ValueError, match="keyword"):
        ParameterDrive("", (Harmonic(1, 1.),))
    with pytest.raises(ValueError, match="at least one"):
        ParameterDrive("phi", ())
    with pytest.raises(ValueError, match="repeated"):
        ParameterDrive("phi", (Harmonic(1, 1.), Harmonic(1, 2.)))
    with pytest.raises(ValueError, match="index"):
        ParameterDrive.cosine("phases", 1., index=-1)
    base = junction()
    build = PeriodicHamiltonian.from_parameters
    with pytest.raises(ValueError, match="ParameterDrive"):
        build(reference_model, base, ["phi"], 1.)
    with pytest.raises(ValueError, match="only once"):
        build(reference_model, base, [ParameterDrive.cosine("phi", .1)]*2, 1.)
    with pytest.raises(ValueError, match="undriven value"):
        build(reference_model, base, ParameterDrive.cosine("rho_ws", .1), 1.)
    with pytest.raises(ValueError, match="real scalar"):
        build(make_model, dict(phases=[0.]), ParameterDrive.cosine("phases", .1), 1.)
    with pytest.raises(ValueError, match="does not select"):
        build(reference_model, base, ParameterDrive.cosine("phi", .1, index=0), 1.)
    with pytest.raises(ValueError, match="construction failed at t="):
        build(reference_model, base, ParameterDrive.cosine("gamma", 1.), 1.)
    with pytest.raises(TypeError, match="return a Hamiltonian"):
        PeriodicHamiltonian.from_callable(lambda t: None, 1.)
    for omega in (0., -1., np.inf, 1j):
        with pytest.raises(ValueError, match="omega"):
            build(reference_model, base, ParameterDrive.cosine("phi", .1), omega)
    with pytest.raises(ValueError, match="three samples"):
        build(reference_model, base, ParameterDrive.cosine("phi", .1), 1., samples=2)
    with pytest.raises(ValueError, match="sampling_tolerance"):
        build(reference_model, base, ParameterDrive.cosine("phi", .1), 1., sampling_tolerance=-1.)
    assert_allclose(fold_quasienergies([-.6, .4, .5, 1.7], 1., center=.2), [.4, .4, .5, -.3], atol=1e-15)
    assert_allclose(quasienergy_distance(.05, .95, 1.), .1, atol=1e-15)


def test_phase_drive_harmonics_follow_jacobi_anger():
    base = junction()
    amplitude, theta, omega = 1.3, .4, .7
    periodic = PeriodicHamiltonian.from_parameters(reference_model, base,
                                                   ParameterDrive.cosine("phi", amplitude, theta), omega, samples=32)
    diagnostics = periodic.metadata["floquet"]
    assert diagnostics["interpolation_error"] < 1e-12
    assert diagnostics["hermiticity_error"] < 1e-15
    assert periodic.metadata["drives"][0]["parameter"] == "phi"
    # The tunnelling factors exp(-+i phi/4) give J_k(A/4) harmonics; H_0 renormalizes Gamma by J_0^2.
    renormalized = reference_model(**{**base, "gamma": base["gamma"]*jv(0, amplitude/4)**2})
    assert difference(periodic.harmonic(0), renormalized.operator) < 1e-14
    static = reference_model(**base)
    tunnelling = static.operator - reference_model(**{**base, "gamma": 0.}).operator
    for k in range(1, 7):
        factor = jv(k, amplitude/4)*np.exp(1j*k*theta)
        if k % 2:
            expected = static.observables["phase_derivative"]*complex(4*(-1)**((k-1)//2)*factor)
        else:
            expected = tunnelling*complex((-1)**(k//2)*factor)
        assert difference(periodic.harmonic(k), expected) < 1e-14
        assert difference(periodic.harmonic(-k), expected.dagger()) < 1e-14
    t = .123
    exact = reference_model(**{**base, "phi": base["phi"] + amplitude*np.cos(omega*t+theta)})
    assert difference(periodic.operator_at(t), exact.operator) < 1e-12
    assert periodic.operator_at(t).hermiticity_error() < 1e-15
    assert periodic.period == pytest.approx(2*np.pi/omega)
    with pytest.raises(ValueError, match="integer"):
        periodic.harmonic(1.)
    # Observables are sampled too: the phase derivative oscillates with the drive.
    current = periodic.observables["phase_derivative"]
    assert set(current) >= {-1, 0, 1}
    assert difference(periodic.time_average().observables["impurity_charge"],
                      static.observables["impurity_charge"]) < 1e-15


def test_linear_parameter_drives_are_band_limited():
    base = junction()
    amplitude, theta, omega = .3, .25, .8
    drives = {name: ParameterDrive.cosine(name, amplitude, theta) for name in ("detuning", "field", "u")}
    results = {name: PeriodicHamiltonian.from_parameters(reference_model, base, drive, omega, samples=12)
               for name, drive in drives.items()}
    for name, periodic in results.items():
        assert periodic.metadata["floquet"]["retained_harmonics"] == [0, 1], name
        assert difference(periodic.harmonic(0), reference_model(**base).operator) < 1e-14
    # H_d contains detuning*(n-1): the drive carries the scalar harmonic -(A/2)exp(i theta).
    scalar = results["detuning"].metadata["floquet"]["scalar_harmonics"]["1"]
    assert_allclose(complex(*scalar), -amplitude/2*np.exp(1j*theta), atol=1e-15)
    expected = (number(0)+number(1)-1)*complex(amplitude/2*np.exp(1j*theta))
    assert difference(results["detuning"].harmonic(1), expected) < 1e-15
    phase = ParameterDrive.cosine("phi", .9)
    combined = PeriodicHamiltonian.from_parameters(reference_model, base, (phase, drives["detuning"]), omega)
    alone = PeriodicHamiltonian.from_parameters(reference_model, base, phase, omega)
    for k in range(1, 6):
        assert difference(combined.harmonic(k), alone.harmonic(k) + results["detuning"].harmonic(k)) < 1e-14


def _delta_model(delta=1.):
    return reference_model(cosh_grid(pairs=1, delta=delta, bandwidth=3.), symmetry=False, compress=False)


@pytest.mark.parametrize("builder, base, drive", [
    (_delta_model, dict(delta=1.), ParameterDrive.cosine("delta", .1)),
    (reference_model, junction(detuning=0., symmetry=True, compress=True), ParameterDrive.cosine("phi", .2)),
    (reference_model, junction(detuning=0., symmetry=True), ParameterDrive.cosine("detuning", .2)),
])
def test_coordinate_guards(builder, base, drive):
    with pytest.raises(ValueError, match="canonical coordinates"):
        PeriodicHamiltonian.from_parameters(builder, base, drive, 1., samples=8)


def test_coordinate_declarations_and_accepted_eta_drive():
    def bare(t):
        return Hamiltonian((1+.1*np.cos(t))*number(0), 1, (1,))
    with pytest.raises(ValueError, match="declare"):
        PeriodicHamiltonian.from_callable(bare, 1.)
    accepted = PeriodicHamiltonian.from_callable(bare, 1., assume_fixed_coordinates=True)
    assert "coordinates" not in accepted.reference.metadata
    # U does not move eta or compressed coordinates.
    base = junction(detuning=0., symmetry=True, compress=True)
    periodic = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("u", .2), 1.)
    assert periodic.eta_labels is not None
    assert periodic.coordinate_fingerprint() == reference_model(**base).coordinate_fingerprint()


def test_sector_and_option_validation():
    def spin_flip(t):
        op = .3*(number(0)-number(1)) + hermitian_pair(.2*np.cos(t)*create(0)*annihilate(1))
        return Hamiltonian(op, 2, (1, -1))
    flipped = PeriodicHamiltonian.from_callable(spin_flip, 1., assume_fixed_coordinates=True)
    with pytest.raises(ValueError, match="S_z"):
        solve_floquet(flipped, 0, Sector(1, 1))
    result = solve_floquet(flipped, 0, Sector(1, None), dict(states=2))
    assert result.quasienergies.shape == (2,)
    with pytest.raises(ValueError, match="parity"):
        PeriodicHamiltonian(1., {0: number(0), 1: create(0)}, 1, (1,))
    with pytest.raises(ValueError, match="k >= 0"):
        PeriodicHamiltonian(1., {0: number(0), -1: number(0)}, 1, (1,))
    with pytest.raises(ValueError, match="FermionOperators"):
        PeriodicHamiltonian(1., {0: number(0), 1: 1.}, 1, (1,))
    with pytest.raises(ValueError, match="undeclared"):
        PeriodicHamiltonian(1., {0: number(0), 1: number(3)}, 1, (1,))
    with pytest.raises(ValueError, match="Hermitian"):
        PeriodicHamiltonian(1., {0: create(0)*annihilate(1)}, 2, (1, -1))
    with pytest.raises(ValueError, match="omega"):
        PeriodicHamiltonian(0., {0: number(0)}, 1, (1,))
    with pytest.raises(ValueError, match="observable"):
        PeriodicHamiltonian(1., {0: number(0)}, 1, (1,), observables={"n": {.5: number(0)}})
    with pytest.raises(ValueError, match="reference must"):
        PeriodicHamiltonian(1., {0: number(0)}, 1, (1,), reference=number(0))
    with pytest.raises(ValueError, match="mode labels"):
        PeriodicHamiltonian(1., {0: number(0)}, 1, (1,), reference=Hamiltonian(number(0), 1, (-1,)))
    assert PeriodicHamiltonian(1., {0: number(0), 2: FermionOperator()}, 1, (1,)).max_harmonic == 0
    with pytest.raises(ValueError, match="PeriodicHamiltonian"):
        solve_floquet(reference_model(**junction()))
    periodic = circular_two_level()
    with pytest.raises(ValueError, match="options"):
        solve_floquet(periodic, 0, Sector(1, None), options=3)
    with pytest.raises(ValueError, match="basis_states"):
        solve_floquet(periodic, 0, Sector(1, None), dict(basis="eigen", states=4, basis_states=2))
    with pytest.raises(ValueError, match="van Vleck order"):
        periodic.van_vleck(3)
    with pytest.raises(ValueError, match="max_harmonic"):
        periodic.van_vleck(1, max_harmonic=-1)
    result = solve_floquet(periodic, 0, Sector(1, None), dict(states=2))
    with pytest.raises(ValueError, match="save_result"):
        save_result(result, "unused.json")


@pytest.mark.parametrize("changes", [
    dict(floquet_cutoff=-1), dict(floquet_cutoff=True), dict(basis="full"), dict(selection="some"), dict(branch="x"),
    dict(method="lanczos"), dict(states=0), dict(steps=1.5), dict(threads=0), dict(tolerance=0.),
    dict(max_memory_gib=np.inf), dict(target_shift=np.nan), dict(zone_center=np.inf),
    dict(tracking_threshold=1.5), dict(step_tolerance=0.), dict(step_check=1),
    dict(drop_scalar_harmonics="yes"), dict(targets=()), dict(targets=(np.nan,)),
    dict(eigenpairs_per_target=8, max_eigenpairs_per_target=4), dict(step_tolerance=1e-6, step_check=False)])
def test_invalid_floquet_options(changes):
    with pytest.raises(ValueError):
        FloquetOptions(**changes)


def test_undriven_limit():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(),
                                                   ParameterDrive.cosine("phi", 0.), .9, samples=8)
    assert periodic.max_harmonic == 0
    sector = Sector(1, 1)
    static = solve(periodic.reference, 2, sector, SolverOptions(method="dense", eigenpairs=4))
    for function, options in ((solve_floquet, dict(floquet_cutoff=3)), (solve_floquet, dict(basis="eigen", basis_states=9)),
                              (floquet_propagator, dict(floquet_cutoff=3, steps=8))):
        result = function(periodic, 2, sector, options)
        assert_allclose(result.quasienergies, static.energies, rtol=0, atol=1e-12)
        assert_allclose(result.mean_energies, static.energies, rtol=0, atol=1e-12)
        center = int(np.flatnonzero(result.harmonic_indices == 0)[0])
        assert_allclose(result.harmonic_weights[:, center], 1., atol=1e-12)
        # Degenerate static states may be rotated; each Floquet state stays in its multiplet.
        multiplet = abs(static.energies[:, None] - static.energies[None, :]) < 1e-9
        assert_allclose(np.sum(result.static_weights[:, :4]*multiplet, axis=1), 1., atol=1e-12)
        assert result.labels.tolist() == [0, 1, 2, 3]


@pytest.mark.parametrize("drop", [True, False])
def test_decoupled_dot_bessel_weights_and_scalar_gauge(drop):
    amplitude, omega = .77, 1.1
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(gamma=0., detuning=.3),
                                                   ParameterDrive.cosine("detuning", amplitude), omega, samples=8)
    for sector in (Sector(0, 0), Sector(1, 1)):
        options = dict(floquet_cutoff=10, states=4, drop_scalar_harmonics=drop)
        sambe = solve_floquet(periodic, 2, sector, options)
        # H(t) commutes at all times; the scheme reduces to Gauss quadrature of the phase.
        propagator = floquet_propagator(periodic, 2, sector, {**options, "steps": 128, "step_check": False})
        static = sambe.static
        assert_allclose(sambe.quasienergies, static.energies[:4], rtol=0, atol=1e-12)
        assert_allclose(sambe.mean_energies, static.energies[:4], rtol=0, atol=1e-12)
        assert_allclose(propagator.quasienergies, static.energies[:4], rtol=0, atol=1e-10)
        # Kept constants give charge n_d-1 in the drive; dropping them leaves n_d.
        charge = np.rint(static.observables["impurity_charge"][sambe.labels])
        charge = charge if drop else charge - 1
        expected = jv(sambe.harmonic_indices[None, :], charge[:, None]*amplitude/omega)**2
        assert_allclose(sambe.harmonic_weights, expected, rtol=0, atol=1e-10)
        expected = jv(propagator.harmonic_indices[None, :], charge[:, None]*amplitude/omega)**2
        assert_allclose(propagator.harmonic_weights, expected, rtol=0, atol=1e-8)


def test_circular_two_level_exact():
    d0, amplitude, omega = 1., .3, 1.2
    periodic = circular_two_level(d0, amplitude, omega)
    assert periodic.metadata["floquet"]["retained_harmonics"] == [0, 1]
    sector = Sector(1, None)
    width = np.hypot(d0-omega, amplitude)
    exact = [(omega+width)/2, (omega-width)/2]
    sambe = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=6, states=2))
    assert match(sambe.quasienergies, exact, omega) < 1e-12
    propagator = floquet_propagator(periodic, 0, sector, dict(floquet_cutoff=6, states=2, steps=64))
    assert match(propagator.quasienergies, exact, omega) < 1e-9
    assert max(propagator.metadata["step_errors"]) < 1e-9
    # Static state 0 is mode 1 occupied (energy -d0/2); E_0+omega is resonant with E_1.
    rwa = rwa_model(periodic, [0, 1], [1, 0], 0, sector)
    assert match(rwa.quasienergies, exact, omega) < 1e-14
    assert rwa.effective_matrix.shape == (2, 2)
    assert rwa.metadata["photons"] == [1, 0]
    sz = number(0) - number(1)
    for order, coefficient in ((0, d0/2), (1, d0/2 - amplitude**2/(4*omega)),
                               (2, d0/2 - amplitude**2/(4*omega) - amplitude**2*d0/(4*omega**2))):
        effective = periodic.van_vleck(order)
        assert difference(effective.operator, coefficient*sz) < 1e-15
        assert effective.metadata["floquet_effective"]["order"] == order


def test_van_vleck_error_scaling():
    sector = Sector(1, 1)
    errors = {order: [] for order in (0, 1, 2)}
    for omega in (40., 80., 160.):
        # Two parameters at one harmonic with a relative phase give a finite first-order term.
        drives = (ParameterDrive.cosine("phase", .8), ParameterDrive.cosine("detuning", .5, .7))
        periodic = PeriodicHamiltonian.from_parameters(lambda **p: single(bandwidth=1.5, **p),
                                                       dict(phase=.5, detuning=.15), drives, omega, samples=16)
        exact = solve_floquet(periodic, None, sector, dict(selection="all", floquet_cutoff=4, method="dense"))
        for order in errors:
            effective = periodic.van_vleck(order)
            assert effective.coordinate_fingerprint() == periodic.coordinate_fingerprint()
            energies = solve(effective, None, sector, SolverOptions(method="dense", eigenpairs=15)).energies
            errors[order].append(match(energies, exact.quasienergies, omega))
    ratios = {order: np.array(values[:-1])/np.array(values[1:]) for order, values in errors.items()}
    assert np.all((ratios[0] > 1.6) & (ratios[0] < 2.5)), ratios
    assert np.all((ratios[1] > 3.) & (ratios[1] < 5.)), ratios
    assert np.all(ratios[2] > 6.), ratios


@pytest.mark.dmrg
def test_van_vleck_effective_hamiltonian_with_dmrg():
    pytest.importorskip("tenpy")
    from qdjj_solver import solve as dispatch
    drives = (ParameterDrive.cosine("phase", .8), ParameterDrive.cosine("detuning", .5, .7))
    periodic = PeriodicHamiltonian.from_parameters(lambda **p: single(bandwidth=1.5, **p),
                                                   dict(phase=.5, detuning=.15), drives, 20., samples=16)
    effective = periodic.van_vleck(2)
    sector = Sector(1, 1)
    exact = dispatch(effective, None, sector, {"eigenpairs": 2, "method": "dense"})
    mps = dispatch(effective, None, sector, {"eigenpairs": 2}, backend="dmrg")
    assert_allclose(mps.energies, exact.energies, rtol=0, atol=2e-10)


def test_sambe_matches_propagator_same_basis():
    omega = .9
    drives = (ParameterDrive.cosine("phi", .6), ParameterDrive.cosine("detuning", .2, .3))
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), drives, omega, samples=32)
    sector = Sector(0, 0)
    sambe = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=8, selection="all", dense_dimension=4000))
    propagator = floquet_propagator(periodic, 2, sector, dict(floquet_cutoff=8, selection="all", steps=128))
    assert len(sambe.quasienergies) == len(propagator.quasienergies) == sambe.basis.dimension
    assert match(sambe.quasienergies, propagator.quasienergies, omega) < 1e-9
    assert max(propagator.metadata["step_errors"]) < 1e-9
    assert_allclose(np.sort(sambe.mean_energies), np.sort(propagator.mean_energies), rtol=0, atol=1e-9)
    tracked = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=8, states=2))
    stepped = floquet_propagator(periodic, 2, sector, dict(floquet_cutoff=8, states=2, steps=128))
    assert_allclose(tracked.quasienergies, stepped.quasienergies, rtol=0, atol=1e-9)
    # The t=0 Floquet mode is the same eigenvector of U(T), and states keep their norm.
    overlaps = abs(np.einsum("dn,dn->n", tracked.mode_at(0.).conj(), stepped.mode_at(0.)))
    assert_allclose(overlaps, 1., atol=1e-8)
    assert_allclose(np.linalg.norm(tracked.state_at(.7), axis=0), 1., atol=1e-10)
    assert_allclose(tracked.fock_modes(), tracked.modes)


def _electron_quasienergies(base, amplitude_phi, amplitude_eps, theta, omega, sector, steps=128):
    """Propagate the exact H(t) in the physical electron space (no Fourier series)."""
    bath, u, gamma = base["bath"], base["u"], base["gamma"]
    hopping, zero = np.sqrt(gamma/(np.pi*bath.rho*2))*np.eye(2), np.zeros((2, 2))
    impurity = Impurity.anderson(u, base["detuning"])
    static = electron_hamiltonian(impurity, [bath, bath], [zero, zero], [0., 0.])

    def lead(l, phase):
        phases = [0., 0.]
        phases[l] = phase
        return electron_hamiltonian(impurity, [bath, bath], [hopping if j == l else zero for j in range(2)],
                                    phases) - static
    modes = 2 + 4*bath.levels
    indices = sector_indices(modes, sector.parity, sector.twice_sz)
    def restrict(matrix):
        return matrix[indices][:, indices].toarray()
    contacts = [restrict((lead(l, 0.) - 1j*lead(l, np.pi))/2) for l in range(2)]
    charge = restrict(evaluate(number(0)+number(1)-1, electron_operators(modes)))
    static = restrict(static)

    def hamiltonian(t):
        phi = base["phi"] + amplitude_phi*np.cos(omega*t)
        result = static + amplitude_eps*np.cos(omega*t+theta)*charge
        for contact, phase in zip(contacts, (-phi/2, phi/2), strict=True):
            term = np.exp(.5j*phase)*contact
            result = result + term + term.conj().T
        return result

    def exponential(matrix, h):
        values, vectors = np.linalg.eigh(matrix)
        return (vectors*np.exp(-1j*h*values)) @ vectors.conj().T
    # Fourth-order Gauss-Legendre Magnus integrator, independent of the solver's scheme.
    period = 2*np.pi/omega
    h = period/steps
    nodes = .5 + np.array([-1., 1.])*np.sqrt(3)/6
    propagator = np.eye(len(indices), dtype=complex)
    for j in range(steps):
        h1, h2 = (hamiltonian((j+c)*h) for c in nodes)
        commutator = h1 @ h2 - h2 @ h1
        generator = (h1+h2)/2 + 1j*np.sqrt(3)/12*h*commutator
        propagator = exponential(generator, h) @ propagator
    return -np.angle(np.linalg.eigvals(propagator))/period


def test_electron_oracle_phase_and_detuning_drive():
    base = junction()
    amplitude_phi, amplitude_eps, theta, omega = .6, .2, .3, .9
    drives = (ParameterDrive.cosine("phi", amplitude_phi), ParameterDrive.cosine("detuning", amplitude_eps, theta))
    periodic = PeriodicHamiltonian.from_parameters(reference_model, base, drives, omega, samples=32)
    sector = Sector(1, 1)
    oracle = _electron_quasienergies(base, amplitude_phi, amplitude_eps, theta, omega, sector)
    propagator = floquet_propagator(periodic, None, sector, dict(selection="all", floquet_cutoff=4, steps=128,
                                                                 step_check=False))
    assert len(oracle) == len(propagator.quasienergies) == 210
    assert match(propagator.quasienergies, oracle, omega) < 1e-8
    sambe = solve_floquet(periodic, None, sector, dict(floquet_cutoff=8, states=3))
    assert sambe.metadata["eigensolver"] == "shift-invert"
    distance = quasienergy_distance(sambe.quasienergies[:, None], oracle[None, :], omega).min(axis=1)
    assert np.all(distance < 1e-8)
    # The leakage bound is rigorous but first order; the actual error is far smaller.
    assert np.all(distance <= sambe.residuals)
    assert np.max(sambe.metadata["window_residuals"]) < 1e-9


def test_eigen_basis_complete_reproduces_fock():
    omega = .9
    drives = (ParameterDrive.cosine("phi", .6), ParameterDrive.cosine("detuning", .2, .3))
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), drives, omega, samples=32)
    sector = Sector(0, 0)
    fock = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=8))
    complete = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=8, basis="eigen", basis_states=42))
    assert complete.representation == "static"
    assert_allclose(complete.quasienergies, fock.quasienergies, rtol=0, atol=1e-10)
    overlaps = abs(np.einsum("mdn,mdn->n", complete.fock_modes().conj(), fock.modes))
    assert_allclose(overlaps, 1., atol=1e-8)
    previous = np.inf
    for basis_states in (10, 20, 26):
        truncated = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=8, basis="eigen", basis_states=basis_states))
        error = abs(truncated.quasienergies - fock.quasienergies)
        assert np.all(error <= truncated.fock_residuals)
        assert np.max(truncated.fock_residuals) < previous
        previous = np.max(truncated.fock_residuals)
    with pytest.raises(ValueError, match="degenerate static multiplet"):
        solve_floquet(periodic, 2, sector, dict(basis="eigen", basis_states=30))


def test_harmonic_truncation_bound():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", 2.),
                                                   .9, samples=32)
    sector = Sector(0, 0)
    reference = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=12, states=3))
    residuals, boundaries = [], []
    for harmonics in (1, 2, 3, 4, 6):
        result = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=harmonics, states=3))
        assert np.all(abs(result.quasienergies - reference.quasienergies) <= result.residuals)
        assert np.all(np.array(result.metadata["window_residuals"]) < 1e-9)
        residuals.append(result.residuals)
        boundaries.append(result.boundary_weights)
    assert np.all(np.diff(residuals, axis=0) < 0)
    assert np.all(np.diff(boundaries, axis=0)[1:] < 0)
    assert np.max(reference.residuals) < 1e-9


def test_gauge_voltage_equivalence():
    """A hopping phase phi(t) equals the level drive epsilon - dphi/dt / 2 for one lead."""
    omega, first, second, theta = .8, .7, .3, .4
    phase = ParameterDrive("phase", (Harmonic(1, first), Harmonic(2, second, theta)))
    periodic = PeriodicHamiltonian.from_parameters(single, dict(phase=.5), phase, omega)
    options = dict(selection="all", steps=256, step_check=False)
    driven = floquet_propagator(periodic, None, Sector(1, 1), options)
    separations = []
    for sign in (1, -1):
        level = ParameterDrive("detuning", (Harmonic(1, sign*omega*first/2, -np.pi/2),
                                            Harmonic(2, sign*omega*second, theta-np.pi/2)))
        gauge = PeriodicHamiltonian.from_parameters(single, dict(detuning=.15), level, omega)
        separations.append(match(driven.quasienergies, floquet_propagator(gauge, None, Sector(1, 1), options)
                                 .quasienergies, omega))
    assert separations[0] < 1e-9
    assert separations[1] > 1e-3


def test_rwa_resonance_splitting():
    base = junction()
    sector = Sector(0, 0)
    reference = reference_model(**base)
    static = solve(reference, 2, sector, SolverOptions(method="dense", eigenpairs=3))
    current = ProjectedOperator(reference.observables["phase_derivative"], static.basis)
    element = abs(current.expectation(static.vectors[:, 0], static.vectors[:, 1]))
    omega = static.energies[1] - static.energies[0]
    gaps = []
    for amplitude in (.02, .01):
        periodic = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", amplitude),
                                                       omega, samples=16)
        rwa = rwa_model(periodic, [0, 1], [0, -1], 2, sector, static=static)
        sambe = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=3, states=2), static=static)
        expected = amplitude*element
        rwa_splitting = quasienergy_distance(*rwa.quasienergies, omega)
        splitting = quasienergy_distance(*sambe.quasienergies, omega)
        assert abs(rwa_splitting/expected - 1) < amplitude**2
        gaps.append(abs(splitting - rwa_splitting)/splitting)
    assert gaps[0]/gaps[1] > 3
    assert np.all(abs(sambe.static_weights[:, :2] - .5) < .05)
    # Resonant targets are omega apart: shift-invert finds replicas that must be merged.
    iterative = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=3, states=2, method="shift-invert",
                                                        eigenpairs_per_target=4), static=static)
    assert_allclose(iterative.quasienergies, sambe.quasienergies, rtol=0, atol=1e-10)
    assert sorted(rwa.labels.tolist()) == [0, 1]
    with pytest.raises(ValueError, match="distinct pairs"):
        rwa_model(periodic, [0, 0], [0, 0], 2, sector)
    with pytest.raises(ValueError, match="one nonnegative"):
        rwa_model(periodic, [0], [0, 1], 2, sector)
    with pytest.raises(ValueError, match="different Fock basis"):
        rwa_model(periodic, [0, 1], [0, -1], 1, sector, static=static)
    with pytest.raises(ValueError, match="eigenpairs"):
        solve_floquet(periodic, 2, sector, dict(states=4), static=static)
    other = solve(reference_model(**{**base, "symmetry": True, "compress": True, "detuning": 0.}), 2, sector)
    with pytest.raises(ValueError, match="canonical coordinates"):
        solve_floquet(periodic, 2, sector, dict(states=1), static=other)
    with pytest.raises(ValueError, match="QP Eigenstates"):
        solve_floquet(periodic, 2, sector, static=object())


def test_eta_sector_union():
    base = junction(detuning=0., symmetry=True)
    periodic = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", .8), .9)
    assert periodic.eta_labels is not None
    options = dict(selection="all", floquet_cutoff=2, steps=64, step_check=False)
    full = floquet_propagator(periodic, 2, Sector(0, 0), options)
    parts = np.concatenate([floquet_propagator(periodic, 2, Sector(0, 0, eta), options).quasienergies
                            for eta in (1, -1)])
    assert len(parts) == len(full.quasienergies)
    assert match(parts, full.quasienergies, .9) < 1e-10


def test_periodic_record_roundtrip():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", .5), .7)
    restored = PeriodicHamiltonian.from_record(periodic.record())
    assert restored.fingerprint() == periodic.fingerprint()
    assert restored.coordinate_fingerprint() == periodic.coordinate_fingerprint()
    assert difference(restored.operator_at(.41), periodic.operator_at(.41)) == 0
    assert restored.metadata["floquet"]["samples"] == 64


def test_shift_invert_tracking_skips_replicas():
    # At the S-S* resonance both dressed states carry weight on S and S*. With one eigenpair per target,
    # the S* target first finds the replica (eps+omega, u_{m-1}) of the state found for S. A per-target
    # weight test accepts it; the search must continue until two distinct Floquet states are assigned.
    model = junction(bath=cosh_grid(pairs=1, delta=1., bandwidth=10.), u=1., gamma=.5, detuning=.2)
    sector = Sector(0, 0)
    static = solve(reference_model(**model), 2, sector, SolverOptions(method="dense", eigenpairs=10))
    omega = float(static.energies[1] - static.energies[0])
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("phi", .4),
                                                   omega, samples=16)
    options = dict(floquet_cutoff=3, states=2)
    dense = solve_floquet(periodic, 2, sector, dict(options, method="dense"), static=static)
    iterative = solve_floquet(periodic, 2, sector, dict(options, method="shift-invert", eigenpairs_per_target=1),
                              static=static)
    assert max(iterative.metadata["eigenpairs_per_target"]) > 1
    assert_allclose(iterative.quasienergies, dense.quasienergies, rtol=0, atol=1e-9)
    assert_allclose(iterative.static_weights, dense.static_weights, rtol=0, atol=1e-8)
    assert min(dense.static_weights[0, 1], dense.static_weights[1, 0]) > .25    # strongly mixed pair


def test_shift_invert_matches_dense_and_guards():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", 1.2),
                                                   .9, samples=32)
    sector = Sector(0, 0)
    dense = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, states=3, method="dense"))
    iterative = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, states=3, method="shift-invert"))
    assert iterative.metadata["eigensolver"] == "shift-invert"
    assert_allclose(iterative.quasienergies, dense.quasienergies, rtol=0, atol=1e-9)
    assert_allclose(iterative.static_weights, dense.static_weights, rtol=0, atol=1e-8)
    targeted = solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, selection="all", method="shift-invert",
                                                       targets=(-.1, .47)))
    assert targeted.metadata["targets"] == [-.1, .47]
    # selection="all" reports the centroid branch, so compare modulo omega.
    distance = quasienergy_distance(targeted.quasienergies[:, None], dense.quasienergies[None, :2], .9)
    assert np.min(distance, axis=0).max() < 1e-9
    with pytest.raises(RuntimeError, match="tracking_threshold"):
        solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, states=3, method="shift-invert",
                                                eigenpairs_per_target=2, max_eigenpairs_per_target=8,
                                                tracking_threshold=.999))
    with pytest.raises(MemoryError):
        solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, method="shift-invert", max_memory_gib=1e-6))
    with pytest.raises(MemoryError):
        solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, method="dense", max_memory_gib=1e-6))
    with pytest.raises(MemoryError, match="max_dimension"):
        solve_floquet(periodic, 2, sector, dict(floquet_cutoff=5, max_dimension=100))
    with pytest.raises(MemoryError, match="dense_dimension"):
        floquet_propagator(periodic, 2, sector, dict(dense_dimension=10))
    with pytest.raises(RuntimeError, match="step_tolerance"):
        floquet_propagator(periodic, 2, sector, dict(floquet_cutoff=1, steps=4, step_tolerance=1e-14))
    with pytest.raises(ValueError, match="steps must exceed"):
        floquet_propagator(periodic, 2, sector, dict(floquet_cutoff=4, steps=8))


def test_array_and_callable_drive_variants():
    omega = .8
    scalar = PeriodicHamiltonian.from_parameters(single, dict(phase=.5), ParameterDrive.cosine("phase", .9), omega)
    bath = cosh_grid(pairs=1, bandwidth=2.5)
    hopping = np.sqrt(.35/(np.pi*bath.rho))*np.eye(2)

    def array_model(phases):
        return make_model(Impurity.anderson(1.3, .15), [bath], [hopping], phases=phases)
    indexed = PeriodicHamiltonian.from_parameters(array_model, dict(phases=[.5]),
                                                  ParameterDrive.cosine("phases", .9, index=0), omega)
    for k in range(scalar.max_harmonic+1):
        assert difference(indexed.harmonic(k), scalar.harmonic(k)) < 1e-15
    base = junction()

    def junction_at(t):
        return reference_model(**{**base, "phi": base["phi"] + .4*np.cos(omega*t)})
    sampled = PeriodicHamiltonian.from_callable(junction_at, omega)
    assert sampled.reference.metadata["floquet"]["kind"] == "time-average"
    assert sampled.coordinate_fingerprint() == reference_model(**base).coordinate_fingerprint()
    with pytest.raises(ValueError, match="increase samples"):
        PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", 20.), omega, samples=4)

    def changing(t):
        h = reference_model(**base)
        if t > 0:
            h.observables.pop("impurity_charge")
        return h
    with pytest.raises(ValueError, match="different observables"):
        PeriodicHamiltonian.from_callable(changing, omega, samples=4)


def test_solver_fallbacks_agree():
    periodic = PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", .7), .9)
    sector = Sector(1, 1)
    options = FloquetOptions(floquet_cutoff=3, states=2, steps=64, step_check=False)
    stored = floquet_propagator(periodic, 2, sector, options)
    # A small memory limit forces the propagator to accumulate harmonics without stored samples.
    streamed = floquet_propagator(periodic, 2, sector, dict(floquet_cutoff=3, states=2, steps=64, step_check=False,
                                                            max_memory_gib=1e-3))
    assert_allclose(streamed.quasienergies, stored.quasienergies, rtol=0, atol=1e-13)
    assert_allclose(streamed.modes, stored.modes, rtol=0, atol=1e-12)
    two_level = circular_two_level()
    tiny = solve_floquet(two_level, 0, Sector(1, None), dict(floquet_cutoff=1, states=2, method="shift-invert"))
    assert tiny.metadata["eigensolver"] == "dense"
    with pytest.raises(ValueError, match="exceeds the basis dimension"):
        rwa_model(periodic, [0, 500], [0, 1], 2, sector)


def test_broad_sidebands_and_scalar_gauge_convergence():
    """Strong gate drive of a proximitized dot: the kept constants centre the sidebands."""
    def model(detuning):
        impurity = Impurity.anderson(1., detuning)
        return Hamiltonian(impurity.operator + hermitian_pair(-.05*create(0)*create(1)), 2, (1, -1))
    periodic = PeriodicHamiltonian.from_callable(lambda t: model(1. + 5.8*np.cos(t)), 1., samples=16,
                                                 assume_fixed_coordinates=True)
    sector = Sector(0, 0)
    kept = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=16, states=2))
    converged = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=26, states=2, drop_scalar_harmonics=True))
    assert_allclose(kept.quasienergies, converged.quasienergies, rtol=0, atol=1e-10)
    # Dropping the constants gives |2> charge 2: at M=16 both states are still found,
    # but the residual bound flags the truncation.
    truncated = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=16, states=2, drop_scalar_harmonics=True))
    error = abs(truncated.quasienergies - converged.quasienergies)
    assert np.all(error <= truncated.residuals) and np.max(truncated.residuals) > 1e-2


def test_resonant_chain_cut_by_a_small_window():
    """At resonance the replicas form a chain; a small window must not duplicate a tracked state."""
    periodic = circular_two_level(1., .05, 1.)
    sector = Sector(1, None)
    width = np.hypot(0., .05)
    exact = [(1.+width)/2, (1.-width)/2]
    for harmonics in (1, 2, 6):
        try:
            result = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=harmonics, states=2))
        except RuntimeError as error:
            assert "replicas" in str(error) or "candidates" in str(error)
            continue
        assert match(result.quasienergies, exact, 1.) <= np.max(result.residuals) + 1e-12
    converged = solve_floquet(periodic, 0, sector, dict(floquet_cutoff=6, states=2))
    assert match(converged.quasienergies, exact, 1.) < 1e-12
