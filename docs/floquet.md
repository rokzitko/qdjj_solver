# Floquet calculations

Periodic modulation of model parameters, such as the impurity level, the
superconducting phase, the hybridization or the interaction, dresses the subgap
states and splits them at resonances. This page explains how to specify a drive,
which method to use, how to read the results, how to compute probe matrix elements
and spectra, and how to check convergence. The calculation yields quasienergies,
harmonic weights, Floquet modes and transition matrix elements. Specify Floquet-state
populations explicitly; predicting steady-state populations requires an additional model of dissipation.

The theory is developed in the physics manual: chapters “Floquet theory for periodically
driven systems” (`ch:floquet`), “Driving the junction” (`ch:driving`), “Numerical methods for driven
problems” (`ch:drivenmethods`) and “Probing the driven junction” (`ch:probes`). The
[tutorials](../examples/README.md) work through examples from an isolated dot to
spectroscopy of the driven junction.

Sambe diagonalization, the one-period propagator and the rotating-wave model use
the QP backend. The high-frequency effective Hamiltonian is an ordinary static
`Hamiltonian`, so it can also be solved with DMRG.

## Conventions

For a first driven calculation, read [Specifying a drive](#specifying-a-drive)
and [Methods](#methods), then use the conventions below to interpret the output.
Start from a checked [static calculation](static_calculations.md), retain fixed
coordinates, and identify the relevant undriven transition before choosing a
frequency. [Probe matrix elements and spectra](spectroscopy.md) supplies a small,
self-contained spectroscopy example.

Energies, the drive frequency $\Omega$ and tolerances are in the units of the
Hamiltonian (typically $\Delta=1$), with $\hbar=1$. The period is
$T=2\pi/\Omega$ and the time origin is $t=0$. A periodic Hamiltonian is stored as
harmonics

```math
H(t)=\sum_k H_k\,e^{ik\Omega t},\qquad H_{-k}=H_k^\dagger,\qquad
H_k=\frac{1}{T}\int_0^T H(t)\,e^{-ik\Omega t}\,dt .
```

Only $k\ge 0$ is stored. A Floquet state and its periodic modes $u_m$ are

```math
\psi(t)=e^{-i\varepsilon t}\sum_m e^{im\Omega t}\,u_m ,
```

and the truncated Sambe (extended-space) eigenvalue problem is

```math
\sum_{m'}\left(H_{m-m'}+m\Omega\,\delta_{mm'}\right)u_{m'}=\varepsilon\,u_m ,
\qquad -M\le m,m'\le M .
```

The Floquet cutoff $M$ is `FloquetOptions.floquet_cutoff`. The pair
$`(\varepsilon+\ell\Omega,\,u_{m-\ell})`$ describes the same physical state for
every integer $\ell$, and moves the harmonic centroid
$\langle m\rangle=\sum_m m\lVert u_m\rVert^2$ to $\langle m\rangle+\ell$. The
mean energy

```math
\bar E=\varepsilon-\Omega\langle m\rangle=\frac{1}{T}\int_0^T\langle\psi(t)\vert H(t)\vert\psi(t)\rangle\,dt
```

is invariant under this choice and is always reported as `mean_energies`.

Quasienergies are reported on one branch:

- `branch="label"`, the default for tracked states: let $a$ be the static state
  assigned to a Floquet state. The branch is chosen so that the centroid of the
  projection onto $a$,
  $c_a=\sum_m m\lvert\langle a\vert u_m\rangle\rvert^2/\sum_m\lvert\langle a\vert u_m\rangle\rvert^2$,
  lies in the interval from $-1/2$ to $1/2$. The quasienergy then tends to the
  static energy as the drive vanishes. It also varies continuously through a
  one-photon resonance and under strong diagonal driving.
- `branch="centroid"`, the default for `selection="all"`: the total centroid
  $\langle m\rangle$ lies in the interval from $-1/2$ to $1/2$.
- `folded` maps every quasienergy into the zone from $c-\Omega/2$ to $c+\Omega/2$,
  where $c$ is `zone_center`.

The branch shift can move mode components beyond $\lvert m\rvert=M$. The stored
harmonic window then grows accordingly, and the `harmonic_indices` field lists the
retained $m$.

The **energy zero** is that of the model ([models](models.md#units-and-energy-reference)).
$H_0$ contains all c-number constants, so quasienergies are directly comparable
with static energies. Driving $\epsilon$ or $U$ also drives the constants of
$H_d=(U/2)(n_d-1)^2+\epsilon(n_d-1)$. A constant term in $H_k$ with $k\ne0$
multiplies $\psi(t)$ by a periodic global phase. Such terms are kept by default,
so the Floquet problem uses $H(t)$ exactly as built; their values are recorded in
`periodic.metadata["floquet"]["scalar_harmonics"]`. With
`drop_scalar_harmonics=True` they are removed. This leaves the quasienergies, mean
energies and static weights unchanged; only the harmonic weights of $u_m$ change.
For a decoupled dot driven by $\epsilon(t)=\epsilon_0+A\cos\Omega t$, a state of
charge $n_d$ has harmonic weights $J_m((n_d-1)A/\Omega)^2$ with the constants kept
and $J_m(n_dA/\Omega)^2$ with them dropped. Keeping them measures the charge from
half filling, which centres the sidebands and needs a smaller Floquet cutoff $M$.

## Specifying a drive

`PeriodicHamiltonian.from_parameters` samples a model builder over one period:

```python
import numpy as np
from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model
from qdjj_solver.qp_solver import solve_floquet

base = dict(bath=cosh_grid(pairs=2, bandwidth=5.), u=1.6, gamma=.3, phi=.8, detuning=.2,
            symmetry=False, compress=False)
drives = [ParameterDrive.cosine("phi", .05),                        # phi0 + 0.05 cos(omega t)
          ParameterDrive.cosine("detuning", .02, phase=np.pi/2)]    # epsilon0 + 0.02 cos(omega t + pi/2)
periodic = PeriodicHamiltonian.from_parameters(reference_model, base, drives, omega=.57)
states = solve_floquet(periodic, cutoff=2, sector=Sector(0, 0), options=dict(floquet_cutoff=4, states=2))
print(states.quasienergies, states.static_weights[:, :2])
```

Each `ParameterDrive` displaces one builder keyword,
$p(t)=p_0+\sum_n A_n\cos(n\Omega t+\theta_n)$, with one `Harmonic(order, amplitude, phase)`
per term. Amplitudes carry the parameter's units, so phase amplitudes are in
radians. `index=` drives one element of an array keyword, such as the `phases`
of `make_model`. Several parameters may be driven together, with relative phases
$\theta_n$. `PeriodicHamiltonian.from_callable(function, omega)` samples an
arbitrary `function(t)` that returns a `Hamiltonian`. The static states that label
the dressed states are eigenstates of a reference Hamiltonian: the undriven model
`builder(**base)` for `from_parameters`, and for `from_callable` the `reference=`
argument or, without it, the time average $H_0$. A builder failure, such as
$\Gamma(t)<0$, is reported with the time at which it occurred.

**Sampling.** The model is built at `samples` equally spaced times (default 64).
Each operator coefficient is Fourier transformed, which retains harmonics up to
$k=\lceil N/2\rceil-1$ for $N$ samples. With `verify=True` the builder is also
evaluated at the midpoints between samples. There the $\ell_1$ norm of the
coefficient differences must not exceed `sampling_tolerance` (default $10^{-9}$).
This norm bounds the operator-norm error caused by aliasing, harmonic truncation
and the removal of small coefficients. `periodic.metadata["floquet"]` records the norms
`harmonic_l1`, `interpolation_error`, `hermiticity_error` and `cleanup_bound`.
Model observables, for example `phase_derivative`, are sampled as well and
stored as `periodic.observables[name][k]` for later Floquet matrix elements.

**Fixed canonical coordinates.** The QP basis is defined in the canonical
(Bogoliubov) coordinates of the model. Every sample must have the same
canonical coordinates as the reference, compared through `coordinate_fingerprint()`,
a checksum of the declared electron-to-canonical map; otherwise construction raises
`ValueError`. The implemented methods assume fixed coordinates; a moving basis
would additionally require the connection term $-iU^\dagger\partial_tU$. The
following drives are therefore rejected:

- the gap $\Delta$ or the bath grid, which define the Bogoliubov transformation;
- a phase drive with paired-mode compression, whose retained modes depend on $\phi$;
- a detuning or $\rho W_N$ drive of `reference_model(symmetry=True)`, which
  switches the eta representation on and off;
- driven phases or direct hopping with `bath_reference="coupled"`.

Build Floquet models with `symmetry=False, compress=False`. The check compares
coordinates rather than options, so a drive that leaves them fixed, such as $U$
in eta coordinates, is accepted. A Hamiltonian without
`metadata["coordinates"]` must be accepted explicitly with
`assume_fixed_coordinates=True`; only its mode labels are then compared.

Each harmonic must also conserve the requested sector, because a projection
silently drops terms that leave it. Drives of $\epsilon$, $\phi$, $\Gamma$, $U$
and a longitudinal field conserve parity and $S_z$, while a transverse field
requires `twice_sz=None`. A detuning drive breaks eta parity.

### Phase drive and voltage

Phases enter the hoppings as $e^{i\phi_l/2}C_l^\dagger T_l d$, so a time-dependent phase
describes a voltage drive. With $\phi=\phi_R-\phi_L$ and the bias voltage $V=V_L-V_R$, the
electric potential of the left lead relative to the right one, the Josephson relation is
$\dot\phi=2eV/\hbar$ ([conventions](reference/conventions.md)). For a weak
drive $\phi(t)=\phi_0+A\cos(\Omega t+\theta)$,
$`H_{\pm1}\approx(A/2)e^{\pm i\theta}\,\partial H/\partial\phi`$: the drive couples through the
`phase_derivative` operator; the current operator is $2e/\hbar$ times it. The exact harmonics
(a Jacobi–Anger series) and the gauge equivalence with a level drive are derived in the
physics manual, section “Phase drive, voltage and the Josephson relation”
(`sec:phase-voltage`).

## Methods

### Sambe diagonalization

`solve_floquet(periodic, cutoff=None, sector=..., options=None, static=None)`
diagonalizes the truncated Sambe matrix.

- `basis="fock"` uses the full QP Fock basis of the sector at the given cutoff,
  so the Sambe dimension is $(2M+1)d$, with $d$ the dimension of the sector.
  Problems up to `dense_dimension` (default
  2000) are diagonalized densely; larger ones use shift-invert.
  - Shift-invert factorizes $Q-\sigma$ by sparse LU at each target $\sigma$ and
    runs ARPACK. The default targets are the static energies of the tracked states.
  - The vectors from all targets are merged by pivoted QR and a Rayleigh–Ritz
    step, and duplicate replicas are removed.
  - Each tracked state must be assigned a distinct Floquet state with at least
    `tracking_threshold` weight on it. Near a resonance, the eigenpairs closest to
    one target can be replicas of the state assigned at another target; the number
    of eigenpairs of a target is then doubled, up to `max_eigenpairs_per_target`.
  - LU fill-in grows quickly with $D$ and $M$. When the estimated memory exceeds
    `max_memory_gib`, the solver raises `MemoryError` and suggests the static basis.
- `basis="eigen"` truncates to the `basis_states` lowest static eigenstates
  of the reference Hamiltonian, uses the matrix elements $\langle a\vert H_k\vert b\rangle$
  between them and is diagonalized densely.
  - Choose the truncation to retain complete degenerate static multiplets;
    splitting a multiplet raises `ValueError`.
  - `fock_residuals` reports each state's residual in the full cutoff Fock
    space, providing a diagnostic of convergence in the number $K$ of retained static states.
- `selection="tracked"` (default) assigns the lowest `states` static states to
  Floquet states by maximizing their static weights. A weight below
  `tracking_threshold` raises `RuntimeError`.
- `selection="all"` returns one replica of every Floquet state, ordered by quasienergy.

`residuals` combine two parts:

- the window residual $\lVert(Q-\varepsilon)u\rVert$, which measures eigensolver
  accuracy and is enforced against `residual_tolerance`;
- the leakage of $H$ into the omitted harmonics $\lvert m\rvert>M$.

Both are also stored separately, as `metadata["window_residuals"]` and
`metadata["leakage"]`. Because the untruncated Sambe operator is
self-adjoint, an exact quasienergy in the chosen many-body space lies within the
combined residual. The bound is first order; for an isolated quasienergy the
actual error is usually much smaller.

A precomputed QP `Eigenstates` result can be passed as `static=` to reuse the
static states in a scan. It is checked for coordinates, Fock basis and count.

### One-period propagator

`floquet_propagator` has the same signature. Like QuTiP's `FloquetBasis`, it
avoids the extended (Sambe) space.

- **Integration.** The Schrödinger equation with $H(t)=\sum_k H_k e^{ik\Omega t}$ is
  integrated over one period, with all stored harmonics and no Floquet cutoff. The
  scheme is the fourth-order commutator-free Magnus integrator, with `steps` steps
  per period.
  Exponentials come from eigendecompositions, so every step is exactly unitary.
- **Quasienergies.** Folded quasienergies are $-\arg\lambda/T$ for the eigenvalues
  $\lambda$ of $U(T)$, obtained by complex Schur decomposition.
- **Modes.** Mean energies from time averages select the central replica. The
  modes are then Fourier transformed and truncated to $\lvert m\rvert\le M$;
  `boundary_weights` reports the weight outside this window.
- **Step check.** With `step_check=True` (default) the propagation is done with
  `steps` and with `2*steps` steps, and the finer result is reported.
  `metadata["step_errors"]` holds the Richardson estimate
  $\lvert\Delta\varepsilon\rvert/15$ of its error, which `step_tolerance` can enforce.

The representation is dense and limited to `dense_dimension`. A propagator in
the Fock basis is an independent check of the Floquet cutoff $M$ for small $D$. In
the static basis it is practical for production calculations.

### High-frequency expansion

`periodic.van_vleck(order)` returns the van Vleck effective Hamiltonian through
order $\Omega^{-n}$, where $n$ is `order` (0, 1 or 2, default 2); the expansion is
given in the physics manual, section “The high-frequency expansion” (`sec:van-vleck`).

The commutators are formed in the full fermionic algebra before any QP
truncation. The result is a static `Hamiltonian` in the same coordinates, which
`solve` can treat with either backend. Its eigenvalues approximate the
quasienergies with error $O(\Omega^{-(n+1)})$ at order $n$. The first-order term
vanishes when $H_{-k}=H_k$, as for a single-parameter cosine drive.

The expansion requires $\Omega$ to exceed the many-body spectral width,
including the bath bandwidth. It provides a check for small models in that regime.
For subgap resonances with $\Omega$ of the order of $\Delta$, use Sambe diagonalization or the propagator. The
number of operator terms grows with products of the harmonics' term counts. The
copied observables are time averages without the $O(1/\Omega)$ kick-operator
correction.

### Rotating-wave model

`rwa_model(periodic, states, photons, ...)` restricts the Sambe matrix to the
pairs $\vert a_i,n_i\rangle$ of static states and photon numbers:

```math
Q^{\mathrm{RWA}}_{ij}=\big(H_{n_i-n_j}\big)_{a_ia_j}+n_i\Omega\,\delta_{ij}.
```

Choose the photon numbers so that the energies $E_{a_i}+n_i\Omega$ are nearly
degenerate. For absorption from $a$ to $b$ with $E_b-E_a\approx\Omega$, use
`states=[a, b]` and `photons=[0, -1]`. At resonance the splitting is
$2\lvert(H_{n_a-n_b})_{ab}\rvert$, which for a weak phase drive is
$A\lvert\langle a\vert\partial_\phi H\vert b\rangle\rvert$. The model neglects
counter-rotating terms (Bloch–Siegert shift) and the AC Stark shifts from all
other states, both of second order in the drive amplitude. `effective_matrix`
holds $Q^{\mathrm{RWA}}$.

At the particle–hole-symmetric point ($\epsilon=0$, $\rho W_N=0$), a phase drive
conserves eta parity, restricting allowed couplings to states within each eta sector. For
the parameters of [`examples/floquet.py`](../examples/floquet.py)
the lowest two even states have opposite eta parity at $\epsilon=0$, so the
example uses $\epsilon=0.2$.

## Transitions and spectra

A weak probe $`X(t)=\sum_kX_k\,e^{ik\Omega t}`$ induces transitions between Floquet states
at the frequencies $\pm(\varepsilon_\alpha-\varepsilon_\beta+q\Omega)$. By the golden rule, the
line strengths are $\lvert X^{(q)}_{\alpha\beta}\rvert^2$, where $X^{(q)}_{\alpha\beta}$ are the
harmonics of the matrix element between the two Floquet states. Positive line frequencies absorb a probe
quantum; negative ones describe stimulated emission. For a complete set of final states the
strengths obey a sum rule. The physics manual derives these results and discusses the
Mollow and Autler–Townes patterns (`ch:probes`).

<!-- docs-test: setup omega_grid = np.linspace(0., 1.5, 1501) -->

```python
states = solve_floquet(periodic, cutoff=2, sector=Sector(0, 0),
                       options=dict(floquet_cutoff=6, selection="all", method="dense"))
table = states.matrix_elements("charge")                    # FloquetTransitions
lines = table.lines(initial=states.dressed_state(0))        # out of the dressed ground state
spectrum = table.spectrum(omega_grid, width=0.002, kind="rate")
odd = solve_floquet(periodic, cutoff=2, sector=Sector(1, 1), options=dict(floquet_cutoff=6, states=2))
tunnelling = odd.matrix_elements("d_up_dagger", states)     # across parity sectors
```

- **Probes.** Probes are given by name (`periodic.probe(name)`), as a `FermionOperator`
  (treated as time independent), or as a dictionary of harmonics. Named observables use
  their sampled harmonics: the current under a phase drive depends on time, while charge,
  spin and electron operators remain time independent in fixed coordinates.
- **Compatibility.** The two results must belong to the same drive (equal checksums
  `PeriodicHamiltonian.fingerprint()`) and the same scalar gauge (`drop_scalar_harmonics`).
  Sectors, QP cutoffs, representations and methods may differ.
- **`FloquetTransitions`.** It holds `photons` ($q$), `elements` $X^{(q)}_{\alpha\beta}$,
  `frequencies`, `strengths`, `total` (the sum-rule value) and `completeness` (the fraction of
  `total` captured by the computed final states); `complete` marks a complete set of final
  states.
- **`lines(initial, weights)`.** It returns `SpectralLines` out of the chosen initial states,
  by default the dressed ground state. Degenerate multiplets are merged, which makes the strengths
  independent of the basis inside a multiplet. Expectation values (same state, $q=0$) are dropped;
  elastic lines at $q\Omega$ are kept.
- **`spectrum(..., kind)`.** It broadens the lines with Lorentzian or Gaussian shapes (`width` is
  the half width at half maximum). The default `kind="signed"` gives $`S(\omega)=\sum_l s_l\,L(\omega-\omega_l)`$, with
  absorption at $\omega>0$ and stimulated emission at $\omega<0$; `kind="rate"` gives $S(\omega)+S(-\omega)$, the golden-rule
  rate at probe frequency $\omega>0$; `kind="net"` gives absorption minus stimulated emission.

Steady-state populations of Floquet states, which determine the initial-state weights, require a
separate model of the dissipative environment. Here the initial states and their weights are
chosen explicitly. A probe of the pumped transition itself, started from one dressed state,
shows the Mollow pattern of absorption, stimulated emission and elastic lines. The
Autler–Townes doublet appears for an initial state outside the pumped pair.

## Results

`FloquetStates` uses one layout for all methods; its fields are listed in the
[result reference](reference/results.md#floquet-states).

`fock_modes()` expresses static-basis modes in the Fock basis, and
`mode_at(t)` and `state_at(t)` return $u(t)$ and $\psi(t)$. Persistence through
`save_result` currently supports static eigenstates. `PeriodicHamiltonian.record()` and
`from_record()` give a JSON-compatible description of the drive.

## Convergence and pitfalls

- **Separate approximations.** Check convergence in the Floquet cutoff $M$
  (leakage, boundary weights, propagator comparison), the static basis $K$
  (`fock_residuals`), the QP cutoff, the bath discretization, the number of
  samples (`interpolation_error`) and the propagator steps (`step_errors`).
  Assess continuum accuracy through bath refinement with the other approximations converged at each size.
- **Discretized continuum.** A finite bath replaces the continuum by discrete
  levels. Replicas of these levels cross the subgap quasienergies at multiphoton
  resonances and produce many narrow avoided crossings. Their positions and gaps
  depend on the discretization. In the continuum limit they represent the
  decay of the dressed subgap state; the closed finite-bath calculation represents
  their effect as a redistribution of spectral weight among discrete states.
  Follow dressed states through their `static_weights` and identify physical features
  by their stability under bath refinement.
- **Degeneracies.** The modes and labels of states inside a degenerate static
  multiplet depend on the basis. Compare whole multiplets. At a resonance both
  dressed states carry weights close to one half on the resonant static states.
- **Resonant chains.** At a resonance the replicas $\vert a,m\rangle$ and $\vert b,m-1\rangle$ form
  a chain through the whole Sambe space. A small Floquet cutoff cuts it and creates edge states that
  perturb the interior replicas at the order of the Rabi coupling. Use several harmonics even
  for weak drives ($M\gtrsim5$). The solver raises an error when two tracked states are
  replicas of one Floquet state.
- **Branches.** Quasienergies that differ by multiples of $\Omega$ describe the
  same state. Branch choices near zone edges and resonances are conventions;
  `folded` and `mean_energies` are the robust quantities. The choice of a lowest
  quasienergy depends on the selected zone; use the label of the dressed static
  ground state to identify that reference state.

## Validation

The tests of the Floquet methods, of probe matrix elements and of spectra, and the
independent electron-space propagation, are listed in
[validation](validation.md#independent-physics-checks).

## References

The Sambe (extended-space) formulation follows Shirley and Sambe; the
high-frequency expansion follows Eckardt and Anisimovas, and Bukov, D'Alessio
and Polkovnikov; the propagator uses the commutator-free Magnus integrators of
Blanes and Moan, and of Alvermann and Fehske. Full citations are listed under
[Floquet methods](references.md#floquet-methods).
