# Numerical methods and convergence

## Separate approximations

Both solvers act on the same finite many-body Hamiltonian. The bath energies
and weights, QP cutoff, MPS bond dimension, number of sweeps, local eigensolver
tolerances, and choice of targeted states each affect the calculation.
Discarding small Hamiltonian coefficients is a further approximation.
`bandwidth` is a physical parameter:
quadrature refinement at fixed bandwidth targets a finite-band continuum model,
not automatically the wide-band limit.

## QP diagonalization

The QP expansion follows the authors' [ASQ variational manuscript](references.md#qp-variational-method),
building on their two-QP treatments of the superconducting Anderson impurity and
quantum-dot Josephson junction. Those method references are distinct from the
software citation; this implementation allows higher QP cutoffs.

The basis retains all impurity configurations and at most `cutoff` occupied
bath modes, subject to parity, spin, eta, and optional particle-number
constraints. Integer cutoffs must satisfy `0 <= cutoff <= h.bath_modes`, counting
the modes remaining after any model compression. `cutoff=None` includes every
configuration of that finite bath; oversized integers are not clamped.
The Hilbert-space dimension is calculated exactly before constructing the basis,
so oversized calculations can be identified before allocating their storage.

Writing $P_q$ for the projector onto the selected symmetry sector and bath
occupations with total QP number at most $q$, the approximation is

```math
H_q=P_q H P_q,\qquad H_q v_i=E_i v_i.
```

The impurity modes are exempt from this occupation cutoff. When comparing
cutoffs, keep the symmetry sector, bath discretization, quasiparticle vacuum,
and canonical basis fixed so that the retained Hilbert spaces are nested.

Dense diagonalization, sparse ARPACK diagonalization, and a matrix-free method
solve $PHP$. The matrix-free method applies the Hamiltonian to vectors without
storing the full matrix. Automatic method selection depends on the Hilbert-space
dimension, with basis and matrix-storage limits checked in advance.
`method='schur'` eliminates a diagonal two-QP block below its poles when seeking
the lowest state at cutoff two.
The implementation treats off-diagonal entries of magnitude at most `1e-13`
in input energy units as numerical roundoff in that block. Larger entries are
rejected. Eliminating the resulting diagonal block is exact, while the accepted
small off-diagonal entries are omitted; the reconstructed state must still pass
the residual check against the complete projected Hamiltonian. Use sparse or
matrix-free diagonalization when those small entries are physically relevant.

Residuals are $\lVert PHPv-Ev\rVert$ inside the projected space. They do not measure
omitted-QP error. Refine both cutoff and bath resolution; nested-cutoff ground
energies obey Rayleigh–Ritz variational monotonicity. `qp_weights[i,q]` is the
probability of finding $q$ bath QPs in normalized eigenstate $i$; these QP-number
components are not individually normalized.

For an observable-resolved study through six QPs, see
[QP-cutoff convergence](qp_convergence.md). It follows energy, local-moment
fraction, impurity charge probabilities, and the total dot–bath spin correlation,
with separate bath and unrestricted-reference checks, archived numerical data,
and a small runnable finite-bath example.

When several eigenstates are calculated iteratively, pivoted QR orthogonalization
and a Rayleigh–Ritz calculation in their span restore orthonormal eigenvectors,
including complex degenerate and near-degenerate multiplets. Rank-deficient
subspaces are rejected. The refined states must still pass the full projected
residual and orthogonality checks.
These checks validate the returned span, not the completeness of a degenerate
multiplet: single-start iterative diagonalization can miss a repeated eigenvalue.
At exact degeneracy, check multiplicities with dense or additional symmetry-resolved
calculations.

## Finite MPS and MPO

The DMRG solver uses TeNPy and groups one or two canonical fermion modes per MPS
site (`group_size`). Conserved quantum numbers follow the mode labels.
`mode_order` specifies their order along the MPS. Construction of the MPO
preserves fermionic operator order and Jordan–Wigner signs, including within
grouped sites. Complex hopping, pairing, and nonlocal interactions are supported.
Ordering modes along an MPS does not transform the bath into a physical chain;
all original couplings remain present.

Exactly proportional MPO channels are combined to reduce the MPO bond
dimension, while preserving conserved quantum numbers and identity operators.
The original channels must be reconstructed exactly in floating-point arithmetic;
there is no singular-value or coefficient threshold. Both the original and
reduced bond dimensions are recorded. If the intermediate matrices would be
too large, this reduction is skipped and the original MPO is used.

DMRG requires `cutoff=None`; a finite QP cutoff is rejected. If the entire model
fits on one MPS site, the requested sector is solved by exact diagonalization.
DMRG includes every mode present in the supplied Hamiltonian. Any dark modes
removed during model construction are already absent from that Hamiltonian.

## Excited states and search shift

Eigenstates are optimized sequentially, orthogonal to the states already found.
Projecting out earlier states can introduce artificial zero eigenvalues in
local problems, which compete with physically positive target energies.

The solver uses the min–max principle to choose an energy shift during the search.
For $k$ requested eigenstates,
the largest Ritz value in any $k$-dimensional trial space bounds the $k$-th exact
eigenvalue from above. Distinct product states form an exactly orthonormal trial
space. Subtracting a value above this bound places all requested eigenvalues below
the artificial zeros, without an extensive shift from empty high-energy modes.
The recorded shift affects optimization only; reported energies and residuals
use the original Hamiltonian.

The shift is also used when seeking a single eigenstate. TeNPy's energy check
divides by $\max(E,1)$, so a negative search energy makes `energy_tolerance` an absolute
per-sweep criterion, independent of extensive empty-mode energies.

Initial product states are chosen to have low local energies and the requested
conserved quantum numbers. By default, each eigenstate is sought from four
initial states: a product state or a previously calculated MPS, followed by
independently randomized states. Local random unitary transformations preserve
the chosen quantum numbers, and their random-number seeds are recorded for
reproducibility. These additional trials can reach conserved subspaces missed
by a particular product state.

A **candidate** is the normalized MPS obtained from one such optimization trial.
For an excited state, its largest absolute overlap with an earlier state,
$\max_j\lvert\langle\psi_j\vert\psi_{\mathrm{candidate}}\rangle\rvert$, must not
exceed `orthogonality_tolerance`. An **excessive-overlap candidate** violates this
test: it retains too much of a state already found, rather than representing a
distinct orthogonal excitation. Such candidates are excluded if any candidate
passes. If all fail, the lowest-energy failed candidate is retained for diagnostics;
"failed" here means failure of this overlap test, not a crashed optimization.
It is not an accepted independent eigenstate. Inspect `metadata["gram_error"]`
and `metadata["finite_problem_converged"]`; `require_convergence=True` raises an
error if the final finite-problem checks fail.

The optional `excitation_operator` specifies a model observable used to generate
an initial fluctuation state $(O-\langle O\rangle)\lvert\psi\rangle$ from the previous
eigenstate. It must preserve the chosen quantum numbers. This can accelerate
the search for a physically motivated excitation; the full variational space
and the other independent trials remain available.

A small residual places the energy close to some eigenvalue; it does not show
that the lowest requested state has been found. Examine trial energies, energy
ordering, overlap matrices, quantum numbers, and calculations with different
initial states or MPS mode orderings. Solve
physically relevant symmetry sectors separately in the same coordinates.

## Worked parity and current example

From the main source directory, run `python -B examples/sector_convergence.py`.
It requires only the base installation. The [complete script](../examples/sector_convergence.py)
constructs a small uncompressed reference junction with `delta=1`, `bandwidth=5`,
`u=1.6`, `gamma=0.3`, and `phi=0.8`, using one mirrored pair of normal-energy nodes
per reservoir.

It solves `Sector(0, None)` and `Sector(1, None)` so **all spin projections** are
included in both parity sectors. It then increases the QP cutoff through
`0, 1, 2, 4, None`. Each sector's lowest energy must decrease or stay constant
within roundoff. Excitation-energy differences and current expectations need
not be monotone under that refinement.

The full finite-space result is approximately:

| Quantity | Value |
|---|---:|
| Lowest even energy / $\Delta$ | `-0.0186381864` |
| Lowest odd energy / $\Delta$ | `-0.2963124911` |
| Ground parity | Odd |
| Parity-changing excitation energy / $\Delta$ | `0.2776743047` |
| Ground-state current / $(2e\Delta/\hbar)$ | `-0.0036057498` |

The gap here is $E_{\text{opposite parity}}-E_{\text{ground parity}}$, evaluated in
the same energy reference. It is a **finite-model parity-changing excitation**, distinct
from a spin excitation within one parity sector. The script checks the current
observable against a centered finite difference of the same sector's energy
under phase changes. This point is away from a parity crossing, so that
derivative is smooth.

The final full-space solve checks the cutoff approximation for this small bath.
For a physical accuracy study, next refine bath nodes at fixed gap, bandwidth,
hybridization, phase, and reservoir measure; converge the solver/cutoff again for
each bath. Stable eigensolver residuals alone do not establish that continuum
limit or justify a continuum binding-energy claim.

## Residuals and convergence checks

For a normalized DMRG state, the full finite-Hamiltonian residual is

```math
r_i=\lVert(H-E_i)\lvert\psi_i\rangle\rVert.
```

It has energy units. For Hermitian $H$ with exact eigenvalues $\lambda_n$,
$\min_n\lvert E_i-\lambda_n\rvert\leq r_i$. This bounds the distance to the spectrum,
but does not by itself bound the error in an observable or identify a unique
eigenvector near degeneracy. A projected QP residual gives the corresponding
statement for $P_qHP_q$ within the retained Hilbert space.

It is computed by expressing $H-E_i$ as an MPO, applying it without compression,
and evaluating the resulting norm by QR factorization within conserved sectors.
No singular values are discarded in this residual calculation.
This avoids subtracting nearly equal $\langle H^2\rangle$ and $\langle H\rangle^2$.
Temporary bonds can be as large as MPS bond dimension times MPO bond dimension.

Local Lanczos convergence is separate: `lanczos_probability_tolerance` sets
TeNPy's dimensionless `P_tol`, the threshold on the estimated probability error
$(r_{\mathrm{Ritz}}/\Delta_{\mathrm{Krylov}})^2$, where $r_{\mathrm{Ritz}}$ is the
local Ritz residual and $\Delta_{\mathrm{Krylov}}$ the estimated local Krylov gap.
For high-accuracy matched-bath comparisons, explicitly set
`lanczos_probability_tolerance=1e-22` with sufficient `chi_max`, and verify
physical residuals; tightening local solves cannot remove bond truncation.

Small discarded weights need not imply small residuals: their square roots set
amplitude errors, which high-energy scales can amplify in the residual.
Vary `chi_max` and local Lanczos precision separately at fixed bath to distinguish
bond truncation from incomplete local optimization. Increasing the bond dimension
can expose a local-solver precision limit; tightening local tolerances cannot
compensate for an insufficient bond dimension. Check full residuals and the
observables of interest for every targeted state. Fixed-bath solver accuracy
does not establish bath convergence.

Saved diagnostics include sweep histories for each eigenstate, actual bond
dimensions, normalization errors, independent trials, residuals, and the full
overlap matrix $G_{ij}=\langle\psi_i\vert\psi_j\rangle$ (the Gram matrix).
`finite_problem_converged` requires stable sweep energies, acceptable residuals,
an overlap matrix close to the identity, and eigenstates ordered by energy.
TeNPy's combined energy/entropy stopping criterion is recorded
separately: eigenvectors can rotate inside a degenerate eigenspace, changing
their entanglement entropy even after their residuals have converged.

`require_convergence=True` stops the calculation with an error if these
finite-problem checks fail. By default, coarse calculations return their
diagnostics so that they can be included in a convergence study.
Without residual evaluation, convergence is unverified. These checks never
establish bath convergence; `bath_convergence_checked` stays false.

## Bond-dimension and bath refinement

For bath generators, parameter selection, and accuracy-versus-time measurements,
see [bath and reservoir representations](bath_representations.md). Its benchmark
results distinguish convergence for a fixed finite bath from convergence toward
the continuum reservoir.
For surrogate fits, `max_nfev` controls function evaluations per start, not
accuracy; see the [evaluation-budget guidance](bath_representations.md#evaluation-budgets-for-high-accuracy-studies)
for why optimizer success and a small held-out kernel error do not certify
continuum observables.

Using `h` and `sector` from the [README calculation](../README.md#first-calculation)
with the `[dmrg]` extra installed:

```python
from qdjj_solver.dmrg_solver import SolverOptions
from qdjj_solver.dmrg_solver.convergence import converge_bond_dimension

result, assessment = converge_bond_dimension(
    h, [32, 64, 128], sector=sector,
    options=SolverOptions(eigenpairs=2, seed_trials=4),
    energy_tolerance=1e-8, observable_tolerance=1e-7,
)
if not assessment["empirical_convergence"]:
    raise RuntimeError("Bond ladder did not establish empirical convergence")
```

Each step starts from the previously calculated MPS eigenstates and includes
independent trials when requested. The helper runs the entire supplied ladder
and returns its **last** result, not the last successful one. It deliberately
sets `require_convergence=False` at each step, even if the supplied options request
strict convergence. Exhausting the ladder without convergence does not itself
raise an error; check `assessment["empirical_convergence"]` as above.
The assessment requires consecutive stable steps (`stable_steps=2` by default).
An optional `callback(result, entry)` can archive each completed step.
`compare_states` reports maximum-overlap assignments, energy/observable changes,
and overlap-matrix singular values. Singular values compare retained subspaces
independently of rotations within degenerate multiplets. Empirical stability is
not converted into an unsupported rigorous continuum error bar.

Then refine the bath with fixed physical parameters, normal-band measure, and
energy reference. Different bath coordinates generally prevent direct MPS
continuation. For gap-edge observables, compare several node resolutions and
possibly an independent positive quadrature supplied through `DiscreteBath`.

Build dissociation thresholds from the lowest allowed many-body energies plus
physical continuum gaps. A discrete bath's lowest level may exceed the physical
edge; using it as the continuum threshold can produce a false binding estimate.
Choose the allowed dissociation channels and required precision for the
particular physical observable being studied.

## Independent validation

Small finite models should agree with unrestricted ED. The numerical checks also
build physical-electron Hamiltonians directly from tensor products, independently
of the solver's Bogoliubov transformation and Fock-basis implementation.
The comparisons include complex contacts and spin–orbit coupling, multiorbital
interactions, permuted modes, positive-energy states, particle-number sectors,
composite observables, and cross-parity transitions.
