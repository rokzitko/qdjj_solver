# Physical models and conventions

## Units and energy reference

Energies, gaps, hybridizations, and energy/residual tolerances use the units
of the supplied Hamiltonian. Both numerical methods use these units directly.
`Impurity.anderson(u, detuning, field)` constructs

```math
H_d=\frac U2(n_d-1)^2+\epsilon(n_d-1)+B S_d^z,
\qquad S_d^z=(n_\uparrow-n_\downarrow)/2.
```

Here $n_d=n_\uparrow+n_\downarrow$. The empty, spin-up, spin-down, and doubly
occupied states have energies $U/2-\epsilon$, $B/2$, $-B/2$, and $U/2+\epsilon$,
respectively, with `detuning` $=\epsilon$ and `field` $=B$.
$B$ is a Zeeman **energy**, supplied as `field`, rather than a magnetic
field in tesla. Positive `field` raises spin-up. For the conventional Anderson
impurity $H_A=\epsilon_d n+U n_\uparrow n_\downarrow+B S_z$, the conversion is

```math
\epsilon=\epsilon_d+\frac U2,\qquad H_d=H_A-\epsilon_d.
```

Thus set `detuning=epsilon_d+u/2`. Add $\epsilon_d$ to the reported impurity energy
to recover the conventional additive reference; excitation differences are
unaffected. For a complete reservoir model, also account for the bath reference
below. Spin observables use units of $\hbar$, and `field` is the energy multiplying
that dimensionless spin.

The functions `make_model` and `reference_model` subtract the disconnected
BCS vacuum energy

```math
E_{\mathrm{BCS}}^{(0)}=
\sum_l\sum_i\left(\xi_{li}-\sqrt{\xi_{li}^2+\Delta_l^2}\right).
```

For one Anderson impurity, the full finite-system energy in the conventional
electron Hamiltonian is therefore

```math
E_{\mathrm{conventional}}=E_{\mathrm{reported}}+\epsilon_d+E_{\mathrm{BCS}}^{(0)}.
```

A coupled-reservoir QP vacuum retains its exact energy difference from that same
reference. Additive constants are retained in the operator expression and
saved Hamiltonian.
A user-built `Hamiltonian` uses exactly its specified additive constant.

## Reservoir measure

The [bath representation guide](bath_representations.md) gives complete parameter
tables, switching examples, convergence properties, and measured accuracy/cost
comparisons for the quadratures and surrogate fits summarized here.

`DiscreteBath(xi, weights, delta, bandwidth)` describes spin-degenerate signed
normal levels. For a flat band of half-width $D$, the normal-state density of
states **per spin** is $\rho=1/(2D)$, and

```math
C_{l\sigma}=\sum_i\sqrt{w_{li}}c_{li\sigma},\qquad w_{li}\simeq\rho_l\,d\xi.
```

Weights represent the normal-state measure, not the BCS density of states.
The contact operator need not be normalized:
$`\{C_{l\sigma},C_{l\sigma}^\dagger\}=\sum_i w_{li}`$. Rescaling the weights to
make this sum one changes the physical tunneling amplitudes and is not done
automatically.
Each signed level contributes two positive-energy Bogoliubov modes. Weights
must be positive. The bath keeps a fixed copy of the supplied energies and
weights; nearly mirrored grids do not qualify for exact symmetry reductions.

`cosh_grid(pairs, delta, bandwidth)` uses Gauss–Legendre quadrature in

```math
\xi=\pm\Delta\sinh t,\quad E=\Delta\cosh t,
\quad 0\le t\le\mathrm{asinh}(D/\Delta).
```

It yields `2*pairs` signed levels **per reservoir** and includes the
$\Delta\cosh t$ Jacobian. Coarse weights are not renormalized; the zeroth-moment
error is `bath.metadata['measure_error']`.

`fit_surrogate(levels, ...)` fits positive residues to the finite-band function

```math
g(i\omega)=\frac{2}{\pi}
\frac{\arctan[D/\sqrt{\omega^2+\Delta^2}]}{\sqrt{\omega^2+\Delta^2}}.
```

Its default logarithmic mesh and absolute least-squares criterion follow the
published construction in [references](references.md). Settings and fit errors
are recorded. A small hybridization-fit error does not certify every excitation.

`chain_expansion(levels, ...)` constructs the Padé ChE bath of Bobok et al.
(2025), with a finite-band target by default and analytic wide-band coefficients
selected by `wide_band=True`. The returned bath has exact mirrored star
coordinates and unnormalized positive contact weights. In wide-band mode,
`bandwidth` only fixes the density normalization and cancels from physical
hybridizations. See [chain-expansion reservoirs](chain_expansion.md) for the
construction, parameters and independent validation.

## Contacts, phases, and spin mixing

The phase convention places the superconducting phases in the hopping terms.
In this gauge, each reservoir has the real-gap Hamiltonian

```math
H_l=\sum_{i,\sigma}\xi_{li}c_{li\sigma}^\dagger c_{li\sigma}
-\Delta_l\sum_i\left(c_{li\uparrow}^\dagger c_{li\downarrow}^\dagger
+c_{li\downarrow}c_{li\uparrow}\right),\qquad \Delta_l>0.
```

This fixes the pairing sign and operator order used for anomalous expectation
values. Its vacuum contribution is included in the subtraction described above.

In `make_model`, `tunneling[l]` is a matrix with two reservoir-spin rows and
`2*impurity.orbitals` impurity-spin-orbital columns. It enters as
$e^{i\phi_l/2}C_l^\dagger T_l d+\mathrm{h.c.}$ Here $C_l$ and $d$ are columns
of reservoir-contact and impurity annihilation operators. `direct[(l,m)]`
specifies $W_{lm}$, for $l<m$, and contributes
$e^{i(\phi_l-\phi_m)/2}C_l^\dagger W_{lm}C_m+\mathrm{h.c.}$ Both matrices may be
complex and spin mixing. Use `twice_sz=None` if spin projection is not conserved.

`bath_reference='coupled'` first diagonalizes the quadratic reservoir Hamiltonian;
that transformation currently requires spin-conserving direct hopping and an
unpolarized quasiparticle vacuum with a nonzero excitation gap. The
isolated-vacuum construction supports general spin-mixing contacts.

`reference_model` constructs a single dot with symmetric lead couplings.
`gamma` is **total** hybridization:
$\Gamma_l=\pi\rho\lvert V_l\rvert^2$, $\Gamma=\sum_l\Gamma_l$.
A symmetric junction uses $\Gamma/2$ per lead and phases $[-\phi/2,+\phi/2]$.
Phases are in radians, with $\phi=\phi_R-\phi_L$. In the `(up, down)` spin order,
the reference model's direct contact matrix is

```math
W_{LR}=W_N I+iW_S\sigma_z
=\begin{pmatrix}W_N+iW_S&0\\0&W_N-iW_S\end{pmatrix},
\qquad W_N=\frac{\texttt{rho\_wn}}{\rho},\quad
W_S=\frac{\texttt{rho\_ws}}{\rho}.
```

Both `rho_ws` and `rho_wn` are dimensionless $\rho W$ inputs. The Hamiltonian term
is $e^{-i\phi/2}C_L^\dagger W_{LR}C_R+\mathrm{h.c.}$ This spin-diagonal $W_S$ term
conserves $S_z$; more general spin-mixing contacts can be supplied to `make_model`.

## Model constructors

Import `reference_model` and `make_model` from `qdjj_solver`. Both return a
`Hamiltonian`, not eigenstates: choose the backend and QP cutoff separately in
`solve`. The name `reference_model` denotes the convenience constructor used for
the reference junction in the authors' [ASQ variational manuscript](references.md#qp-variational-method),
not a separate numerical method. The earlier single-reservoir and QD-junction
variational papers are listed there as well.

### `reference_model`

| Parameter | Default | Meaning and constraints |
|---|---|---|
| `bath` | Required | One `DiscreteBath`, reused for both leads in a symmetric junction |
| `u` | `2.0` | Impurity interaction $U$, in energy units |
| `gamma` | `0.4` | Nonnegative **total** hybridization, in energy units; split equally between junction leads |
| `phi` | `3*pi/5` | Junction phase difference in radians; ignored for a single reservoir |
| `rho_ws` | `0.0` | Dimensionless spin-dependent direct contact $\rho W_S$ defined above |
| `rho_wn` | `0.0` | Dimensionless spin-independent direct contact $\rho W_N$ defined above |
| `geometry` | `"junction"` | `"junction"` for two leads or `"single"` for one; the latter requires both direct contacts to be zero |
| `detuning` | `0.0` | $\epsilon_d+U/2$, in energy units; zero is half filling |
| `field` | `0.0` | Local Zeeman energy; positive values raise spin-up |
| `symmetry` | `True` | Use eta coordinates when eligible; otherwise keep physical-impurity coordinates |
| `compress` | `True` | Freeze dark modes in their vacuum when eligible; use `False` for the full finite-bath spectrum |
| `bath_reference` | `"isolated"` | `"isolated"` BCS reservoirs or `"coupled"` quadratic environment |
| `coefficient_tolerance` | `1e-14` | Nonnegative Hamiltonian-coefficient cleanup tolerance, in energy units |

All scalar physical parameters and the cleanup tolerance must be finite. Eta
coordinates require a junction with an exactly mirrored bath, `detuning=0`,
`rho_wn=0`, and `bath_reference="isolated"`. Compression additionally requires
`rho_ws=0`. Inspect `h.metadata["paired_mode_compression"]` for the actual choice.
These are Python defaults; the unified CLI instead defaults to `compress=false`.

### `make_model`

This constructor supports multiple impurity orbitals and distinct reservoirs.
Let `L` be the number of reservoirs and `M=impurity.orbitals`.

| Parameter | Default | Meaning and constraints |
|---|---|---|
| `impurity` | Required | An `Impurity` with `M` positive-integer orbitals and a Hermitian fermion operator |
| `baths` | Required | Nonempty sequence of `L` `DiscreteBath` objects |
| `tunneling` | Required | `L` finite complex matrices of shape `(2, 2*M)`, in energy units; reservoir-spin rows, impurity-spin-orbital columns |
| `phases` | `None` | `L` finite phases in radians; `None` means all zero |
| `direct` | `None` | Mapping `(l, m)` to finite complex `(2, 2)` contact matrices in energy units, with `0 <= l < m < L`; `None` means no direct contacts |
| `phase_velocities` | `None` | `L` finite phase derivatives with respect to the chosen parameter; `None` means zero, so no `phase_derivative` observable is added |
| `eta_basis` | `False` | Request verified symmetric single-dot eta coordinates; incompatible models raise an error |
| `compress_pairs` | `False` | Freeze dark modes; requires eligible eta coordinates and zero direct hopping |
| `bath_reference` | `"isolated"` | `"isolated"` or `"coupled"`; coupled coordinates require spin-conserving direct hopping and a gapped, unpolarized reservoir vacuum |
| `coefficient_tolerance` | `1e-14` | Finite nonnegative Hamiltonian-coefficient cleanup tolerance, in energy units |
| `direct_derivatives` | `None` | Mapping of observable names to contact mappings in the same form as `direct`; names must not collide with impurity observables or `phase_derivative` |

Eta coordinates require two identical mirrored baths, a half-filled Anderson
impurity, equal real spin-independent tunnel contacts, and only the pure $W_S$
direct contact defined above. They are unavailable with `bath_reference="coupled"`.
See [coordinates and compression](bath_representations.md#coordinates-compression-and-chains)
for representation choices and [input files](interfaces.md#input-files-and-terminal-commands)
for the JSON equivalents.

## General impurities and finite geometries

Impurity modes are `(orbital 0 up, down, orbital 1 up, down, ...)`; their
Hilbert space is untruncated. `Impurity(orbitals, operator)` accepts arbitrary
fermionic interactions, exchange, pair hopping, and local pairing.
`Impurity.from_integrals(one_body, interaction=None, constant=0.0)` uses

```math
H=C+\sum_{ab}h_{ab}d_a^\dagger d_b
+\frac14\sum_{abcd}U_{abcd}d_a^\dagger d_b^\dagger d_d d_c,
```

with $`h_{ab}=\texttt{one\_body}[a,b]`$, antisymmetrized two-electron integrals,
and $C=\texttt{constant}$. All indices label spin-orbitals. The factor $1/4$
accounts for antisymmetry in
both index pairs.

For arbitrary chains, stars, or imported discretizations, use
`Hamiltonian(operator, nimp, spins, eta_labels=None, observables={}, metadata={})`.
`spins[i]` is a twice-spin mode label `+1` or `-1`; the first `nimp` modes are
exempt from the QP cutoff. The operator expression specifies all hoppings
and interactions and hence the physical geometry.

## Fermionic algebra

Monomials use signed, one-based indices: `+(i+1)` creates, `-(i+1)` annihilates.
Products have mathematical order, with the rightmost factor acting first.
The stored normal-ordered form lists creation operators in ascending mode order
and annihilation operators in descending order. Saved terms contain this
operator sequence and the real and imaginary parts of its coefficient.

The Python functions `create(i)`, `annihilate(i)`, and `number(i)` instead take
zero-based mode indices. For example, `number(0).to_records()` returns
`[{"operators": [1, -1], "real": 1.0, "imag": 0.0}]`.

Operator multiplication and changes of basis apply anticommutation relations
to the full operator **before** projection. For a bath-QP annihilation operator
$c$ and a projector $P$ onto its vacuum, $Pc c^\dagger P=P$, whereas
$(PcP)(Pc^\dagger P)=0$. This matters for composite observables at cutoff boundaries.

`Hamiltonian` checks mode ranges, parity, Hermiticity, and requested symmetries.
These symmetry checks also apply to small coefficients. If small coefficients
are discarded when constructing a model, their summed magnitudes are recorded
as an upper bound on the operator-norm change. Conversion to a matrix product
operator (MPO) introduces no further coefficient cutoff.

## Sectors and canonical coordinates

Both solvers use `Sector(parity=1, twice_sz=1, eta=None, particle_number=None)`.
Parity is zero for even, one for odd. `twice_sz=None` retains all spin projections.
`eta` is the product of occupied modes' assigned eta labels and is available only
in compatible representations. `particle_number` fixes total occupation of the
**chosen canonical modes**, when conserved; it is not implicitly physical
electron number in a Bogoliubov basis.

Model constructors save the electron-to-canonical transformation in
`metadata['coordinates']`.
`h.physical_operator(op)` transforms an expression in electron operators using
that fixed map. Electron ordering is impurity modes, followed by reservoir, level, and
spin indices. User-built models may specify the transformation explicitly.
If none is specified, compatibility checks require an identical Hamiltonian
before reusing states in another calculation.

Eta/paired-mode reduction is restricted to verified symmetric single-dot models.
**`compress=True` fixes discarded dark modes in their vacuum.** It may omit
continuum excitations. DMRG retains all modes present in the model, but cannot
restore modes removed by a model constructor. Use `compress=False` for a full finite-bath
spectrum; the `qdjj-solver` terminal command uses this setting for both methods.

Transforming an arbitrary electron operator with `h.physical_operator` requires
an uncompressed model: multiplying separately projected electron operators
would miss vacuum contractions of the discarded modes. The supplied observables
include these contractions. Saved model information identifies the restriction
and the energies of discarded modes.

## Built-in observables

`h.observables` maps names to operators in the model's canonical coordinates.
Solvers report their expectation values as `result.observables[name][i]` for
eigenstate `i`. `Impurity` supplies the default impurity observables; specify
its `observables` dictionary to add operators or replace these definitions.

| Name | Definition | Availability | Units |
|---|---|---|---|
| `impurity_charge` | Total physical impurity electron occupation | Default for every `Impurity` | Electron count; multiply by $-e$ for electric charge |
| `impurity_spin_z` | Sum of $(n_\uparrow-n_\downarrow)/2$ over impurity orbitals | Default for every `Impurity` | $\hbar$ units |
| `double_occupancy_0` | $n_{0\uparrow}n_{0\downarrow}$ | First impurity orbital; subsequent orbitals use `_1`, `_2`, ... | Dimensionless probability |
| `phase_derivative` | $\partial H/\partial\phi$ in a reference junction; $O_x$ defined below for supplied phase velocities in `make_model` | Reference junctions, or nonzero `phase_velocities` | Energy for dimensionless phase; energy per unit $x$ in general |
| `rho_ws_derivative` | $\partial H/\partial(\rho W_S)$ | Reference junctions, including at `rho_ws=0` | Energy |

The `geometry="single"` reference model has the three impurity observables and
no interlead/phase derivative. General models can supply additional impurity
observables and named `direct_derivatives`. Their units follow the supplied
operator coefficients. Built-in operators already include the transformations
and contractions associated with the constructor's representation.

## Physical derivatives and transitions

To vary the lead phases along a parameter $x$, supply
$\mathrm{d}\phi_l/\mathrm{d}x$ as `phase_velocities`. With all other physical
parameters fixed, the resulting `phase_derivative` is

```math
O_x=\sum_l\frac{\mathrm{d}\phi_l}{\mathrm{d}x}
\frac{\partial H}{\partial\phi_l}.
```

This specifies the derivative due to the changing phases; it does not also
differentiate tunneling amplitudes, gaps, or other model parameters.
`direct_derivatives[name][(l,m)]` specifies a physical contact observable even
when the corresponding Hamiltonian coupling is zero. Operators are constructed
in the physical gauge and then transformed. Differentiating a parameter-dependent
basis without connection terms would give a different operator. For phase bias,
multiply $\partial H/\partial\phi$ by $2e/\hbar$ to obtain current.

For an exact eigenstate on a differentiable energy branch, the
Hellmann–Feynman relation gives

```math
I_i=\frac{2e}{\hbar}\left\langle\psi_i\left|
\frac{\partial H}{\partial\phi}\right|\psi_i\right\rangle
=\frac{2e}{\hbar}\frac{\partial E_i}{\partial\phi}.
```

At zero temperature, use the lowest branch for the equilibrium current.
At a level crossing, take the appropriate one-sided branch derivatives.
For approximate states, check agreement between the current expectation and
the energy derivative as part of convergence. A parameter-dependent truncated
variational space can contribute additional terms to the derivative of its
approximate energy.

If the numerical Hamiltonian is expressed in an energy unit $E_0$ (`E0`), the current is
`I = (2*e*E0/hbar) * result.observables["phase_derivative"][i]`, with $e>0$.
Convert `E0` to joules to obtain amperes. With gap `delta=1`, the numerical
expectation is a current in units of $2e\Delta/\hbar$. This defines the current
sign through the stated phase-bias convention.

Transition matrix elements require states and operators in the same canonical
basis. DMRG also checks the MPS mode ordering, site grouping, and conserved
quantum numbers, and includes the required Jordan–Wigner strings for operators
of odd fermion parity. Independently calculated eigenvectors have arbitrary
phases; compare squared matrix-element magnitudes (transition strengths), or
complete degenerate subspaces, when comparing calculations.
