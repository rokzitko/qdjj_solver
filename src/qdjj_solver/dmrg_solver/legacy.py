"""Compatibility implementations for earlier single-dot DMRG entry points.

Both the physical-electron chain and the diagonalized-reservoir QP star are
exact coordinates at fixed discretization. DMRG has its own Schmidt-rank and
sweep convergence parameters, used to establish reference values beyond the
sizes accessible to unrestricted ED. No bath-QP cutoff is imposed here.
"""

from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits


def chain_coefficients(bath):
    """Lanczos tridiagonalization starting at the normalized contact orbital.

    Full reorthogonalization is used. The returned contact norm must multiply
    dot--lead hopping (and both contact norms multiply interlead hopping).
    """
    norm = np.sqrt(bath.weights.sum())
    vector = np.sqrt(bath.weights)/norm
    vectors, diagonal, hopping = [], [], []
    for i in range(bath.levels):
        vectors.append(vector)
        work = bath.xi*vector
        diagonal.append(float(np.dot(vector, work)))
        work -= diagonal[-1]*vector
        if i:
            work -= hopping[-1]*vectors[-2]
        # Two passes protect orthogonality when the normal energies span many decades.
        q = np.column_stack(vectors)
        for _ in range(2):
            work -= q @ (q.T @ work)
        beta = np.linalg.norm(work)
        if i == bath.levels-1 or beta < 2e-14*max(1., np.max(abs(bath.xi))):
            break
        hopping.append(float(beta))
        vector = work/beta
    diagonal, hopping = np.array(diagonal), np.array(hopping)
    matrix = np.diag(diagonal)
    if hopping.size:
        matrix += np.diag(hopping, 1)+np.diag(hopping, -1)
    xi = np.linalg.eigvalsh(matrix)
    vacuum = float(np.sum(xi-np.hypot(xi, bath.delta)))
    return diagonal, hopping, float(norm), vacuum


def dmrg_chain_reference(bath, u=2., gamma=.4, phi=3*np.pi/5, rho_ws=0., geometry="junction",
                   twice_sz=1, chi=128, max_sweeps=40, energy_tolerance=1e-11, threads=1,
                   calculate_variance=True):
    """Unrestricted finite-bath result, using optional TeNPy (``[reference]``).

    This reference routine treats the symmetric single-dot model at half
    filling. The general multi-orbital QP solver does not depend on TeNPy.
    """
    import tenpy
    from tenpy.algorithms import dmrg
    from tenpy.models.model import CouplingMPOModel
    from tenpy.networks.mps import MPS
    from tenpy.networks.site import SpinHalfFermionSite

    if geometry not in ("single", "junction") or twice_sz not in (-1, 0, 1):
        raise ValueError("reference geometry or spin sector not supported")
    if geometry == "single" and rho_ws != 0:
        raise ValueError("interlead hopping requires two reservoirs")
    if chi < 4 or max_sweeps < 4 or energy_tolerance <= 0:
        raise ValueError("invalid DMRG convergence parameters")
    start = perf_counter()
    onsite, hopping, norm, vacuum = chain_coefficients(bath)
    levels = len(onsite)
    leads = 1 if geometry == "single" else 2
    length, dot = leads*levels+1, 0 if leads == 1 else levels
    positions = [list(range(1, levels+1))] if leads == 1 else [list(range(levels-1, -1, -1)), list(range(levels+1, length))]
    phases = [0.] if leads == 1 else [-phi/2, phi/2]
    velocities = [0.] if leads == 1 else [-.5, .5]
    coupling = np.sqrt(gamma/(np.pi*bath.rho*leads))*norm
    spin_hopping = [1j*rho_ws/bath.rho*norm**2, -1j*rho_ws/bath.rho*norm**2]
    derivative_terms = []

    class JunctionChain(CouplingMPOModel):
        def init_sites(self, model_params):
            return SpinHalfFermionSite(cons_N="parity", cons_Sz="Sz")

        def init_terms(self, model_params):
            self.add_onsite_term(u/2, dot, "Id")
            self.add_onsite_term(-u/2, dot, "Ntot")
            self.add_onsite_term(u, dot, "NuNd")
            for l, chain in enumerate(positions):
                for j, site in enumerate(chain):
                    self.add_onsite_term(onsite[j], site, "Ntot")
                    self.add_local_term(-bath.delta, [("Cdu", [site, 0]), ("Cdd", [site, 0])], plus_hc=True)
                    for creation, annihilation in (("Cdu", "Cu"), ("Cdd", "Cd")):
                        if j:
                            self.add_local_term(hopping[j-1], [(creation, [chain[j-1], 0]),
                                                               (annihilation, [site, 0])], plus_hc=True)
                        else:
                            value = coupling*np.exp(.5j*phases[l])
                            self.add_local_term(value, [(creation, [site, 0]), (annihilation, [dot, 0])], plus_hc=True)
                            derivative_terms.append((.5j*velocities[l]*value, [(creation, site), (annihilation, dot)]))
            if leads == 2:
                for spin, (creation, annihilation) in enumerate((("Cdu", "Cu"), ("Cdd", "Cd"))):
                    value = spin_hopping[spin]*np.exp(-.5j*phi)
                    term = [(creation, [positions[0][0], 0]), (annihilation, [positions[1][0], 0])]
                    self.add_local_term(value, term, plus_hc=True)
                    derivative_terms.append((-.5j*value, [(creation, positions[0][0]), (annihilation, positions[1][0])]))

    with threadpool_limits(limits=threads):
        model = JunctionChain(dict(L=length, lattice="Chain", bc_MPS="finite", bc_x="open"))
        sites = model.lat.mps_sites()
        product = []
        for i, site in enumerate(sites):
            if i == dot and twice_sz:
                product.append("up" if twice_sz == 1 else "down")
            else:
                vector = np.zeros(4)
                vector[site.state_index("empty")] = vector[site.state_index("full")] = 1/np.sqrt(2)
                product.append(vector)
        psi = MPS.from_product_state(sites, product, bc="finite", dtype=complex, permute=False,
                                     unit_cell_width=length)
        options = dict(mixer=True, mixer_params=dict(amplitude=1e-5, decay=2., disable_after=8),
                       max_sweeps=max_sweeps, min_sweeps=10, N_sweeps_check=2,
                       max_E_err=energy_tolerance, max_S_err=1e-9,
                       trunc_params=dict(chi_max=chi, svd_min=1e-14, trunc_cut=None),
                       lanczos_params=dict(N_min=3, N_max=60, P_tol=1e-13, E_tol=energy_tolerance/10))
        if length <= 5:
            options.pop("lanczos_params", None)
        engine = dmrg.TwoSiteDMRGEngine(psi, model, options)
        energy, psi = engine.run()
        variance = float(model.H_MPO.variance(psi, exp_val=energy).real) if calculate_variance else None
        spin = float(psi.expectation_value("Sz", [dot])[0].real)
        charge = float(psi.expectation_value("Ntot", [dot])[0].real)
        double = float(psi.expectation_value("NuNd", [dot])[0].real)
        derivative = float(sum(2*(value*psi.expectation_value_term(term)).real for value, term in derivative_terms))
    stats = engine.sweep_stats
    return dict(method="unrestricted-chain-DMRG", tenpy=tenpy.__version__,
                sweep_converged=bool(engine.is_converged()),
                energy=float(energy-leads*vacuum), energy_before_vacuum_subtraction=float(energy),
                bcs_vacuum_energy=leads*vacuum, variance=variance,
                impurity_spin_z=spin, impurity_charge=charge, double_occupancy_0=double,
                phase_derivative=derivative, twice_sz=twice_sz,
                chi_limit=chi, chi_actual=int(max(psi.chi)), sweeps=int(engine.sweeps),
                last_energy_change=float(stats["Delta_E"][-1]),
                last_discarded_weight=float(stats["max_trunc_err"][-1]),
                norm_error=float(np.linalg.norm(psi.norm_test())),
                seconds=perf_counter()-start, bath=bath.record(), u=u, gamma=gamma, phi=phi,
                rho_ws=rho_ws, geometry=geometry, energy_tolerance=energy_tolerance,
                max_sweeps=max_sweeps, threads=threads,
                sweep_energies=[float(e-leads*vacuum) for e in stats["E"]])


def dmrg_reference(bath, u=2., gamma=.4, phi=3*np.pi/5, rho_ws=0., geometry="junction",
                   twice_sz=1, chi=64, max_sweeps=60, energy_tolerance=1e-11, threads=1,
                   calculate_variance=True, calculate_current=True):
    """Unrestricted DMRG in QPs of the exactly diagonalized quadratic bath.

    No QP-number cutoff is made. The bath vacuum is a product in these
    coordinates, which substantially reduces the required Schmidt rank for
    weakly hybridized superconducting impurities. The isolated BCS energy
    reference is retained, including the exact background-junction shift.
    """
    import tenpy
    from tenpy.algorithms import dmrg
    from tenpy.models.model import CouplingMPOModel
    from tenpy.networks.mps import MPS
    from tenpy.networks.site import SpinHalfFermionSite
    from ..common.models import reference_model

    if twice_sz not in (-1, 0, 1) or chi < 4 or max_sweeps < 10:
        raise ValueError("unsupported reference sector or invalid DMRG parameters")
    start = perf_counter()
    with threadpool_limits(limits=threads):
        use_active_modes = rho_ws == 0 and geometry == "junction" and twice_sz != 0
        h = reference_model(bath, u=u, gamma=gamma, phi=phi, rho_ws=rho_ws,
                            geometry=geometry, bath_reference="isolated" if use_active_modes else "coupled")
        length = len(h.spins)//2
        def local_term(key, lattice=True):
            terms = []
            for op in key:
                i = abs(op)-1
                name = ("Cdu" if i%2 == 0 else "Cdd") if op > 0 else ("Cu" if i%2 == 0 else "Cd")
                terms.append((name, [i//2, 0] if lattice else i//2))
            return terms
        class QuasiparticleStar(CouplingMPOModel):
            def init_sites(self, model_params):
                return SpinHalfFermionSite(cons_N="parity", cons_Sz="Sz")

            def init_terms(self, model_params):
                for key, value in h.operator.terms.items():
                    if key:
                        self.add_local_term(value, local_term(key))
                    else:
                        self.add_onsite_term(value, 0, "Id")
        model = QuasiparticleStar(dict(L=length, lattice="Chain", bc_MPS="finite", bc_x="open"))
        product = ["up" if twice_sz == 1 else ("down" if twice_sz == -1 else "empty")]+["empty"]*(length-1)
        psi = MPS.from_product_state(model.lat.mps_sites(), product, bc="finite", dtype=complex,
                                     unit_cell_width=length)
        options = dict(mixer=True, mixer_params=dict(amplitude=1e-4, decay=2., disable_after=12),
                       max_sweeps=max_sweeps, min_sweeps=14, N_sweeps_check=2,
                       max_E_err=energy_tolerance, max_S_err=1e-9,
                       trunc_params=dict(chi_max=chi, svd_min=1e-14, trunc_cut=None),
                       lanczos_params=dict(N_min=3, N_max=80, P_tol=1e-13, E_tol=energy_tolerance/10))
        if length <= 5:
            options.pop("lanczos_params", None)
        engine = dmrg.TwoSiteDMRGEngine(psi, model, options)
        energy, psi = engine.run()
        variance = float(model.H_MPO.variance(psi, exp_val=energy).real) if calculate_variance else None
        derivative = None
        if calculate_current and "phase_derivative" in h.observables:
            derivative = 0j
            for key, value in h.observables["phase_derivative"].terms.items():
                derivative += value*(psi.expectation_value_term(local_term(key, False)) if key else 1.)
            derivative = float(derivative.real)
        def expectation(operator):
            return float(sum(value*(psi.expectation_value_term(local_term(key, False)) if key else 1.)
                             for key, value in operator.terms.items()).real)
        spin = expectation(h.observables["impurity_spin_z"])
        charge = expectation(h.observables["impurity_charge"])
        double = expectation(h.observables["double_occupancy_0"])
    stats = engine.sweep_stats
    return dict(method="unrestricted-bath-QP-DMRG", tenpy=tenpy.__version__,
                energy=float(energy), variance=variance, sweep_converged=bool(engine.is_converged()),
                impurity_spin_z=spin, impurity_charge=charge, double_occupancy_0=double,
                phase_derivative=derivative, twice_sz=twice_sz,
                chi_limit=chi, chi_actual=int(max(psi.chi)), sweeps=int(engine.sweeps),
                last_energy_change=float(stats["Delta_E"][-1]),
                last_discarded_weight=float(stats["max_trunc_err"][-1]),
                norm_error=float(np.linalg.norm(psi.norm_test())),
                seconds=perf_counter()-start, bath=bath.record(), u=u, gamma=gamma, phi=phi,
                rho_ws=rho_ws, geometry=geometry, energy_tolerance=energy_tolerance,
                max_sweeps=max_sweeps, threads=threads,
                bath_vacuum_shift=h.metadata["bath_vacuum_shift"],
                bath_modes=h.bath_modes, paired_mode_compression=h.metadata["paired_mode_compression"],
                energy_reference=h.metadata["energy_reference"],
                sweep_energies=[float(e) for e in stats["E"]])
