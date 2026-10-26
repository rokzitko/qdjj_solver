# Using the solver and saving results

This chapter uses the small installation-test Hamiltonian. Its Python examples
form one session. Begin with the complete setup below; the larger reference
junction is developed in [static calculations](static_calculations.md).

<!-- BEGIN generated: interface-setup (tools/gen_docs.py) -->
```python
from qdjj_solver import Sector, cosh_grid, reference_model, solve


bath = cosh_grid(pairs=1, delta=1.0, bandwidth=5.0)
h = reference_model(bath, u=1.6, gamma=0.3, phi=0.8, compress=False)
sector = Sector(parity=1, twice_sz=1)
qp = solve(h, cutoff=2, sector=sector, backend="qp")
print(qp.energies, qp.qp_weights)

assert abs(qp.energies[0] - (-0.2950595967187083)) < 1e-10
assert max(qp.residuals) < 1e-9
assert abs(qp.qp_weights[0].sum() - 1) < 1e-12
```
<!-- END generated: interface-setup -->

For the DMRG examples, install the `dmrg` extra and add:

<!-- docs-test: requires dmrg -->
```python
mps = solve(h, sector=sector, backend="dmrg",
            options={"eigenpairs": 2, "chi_max": 32, "require_convergence": True})
```

## Calculations in Python

Import model constructors and the functions `solve`, `save_result`, and
`load_result` from `qdjj_solver`. Their arguments are:

```text
solve(hamiltonian, cutoff=None, sector=Sector(), options=None, initial=None, *, backend="qp")
save_result(result, path, save_states=False)
load_result(path)
```

- `solve` returns the calculated `Eigenstates`. Supply the physical `Hamiltonian`
  as `hamiltonian` and the conserved quantum numbers as `sector`. Choose the
  method with `backend="qp"` or `"dmrg"`. `options=None` uses default numerical
  settings; a dictionary changes individual settings. A `SolverOptions` object
  must come from the chosen solver's module.
- `save_result` saves energies, observables, calculation settings, and diagnostics
  to a filename ending in `.json`. Set `save_states=True` to save the wavefunctions
  in a companion file. The function returns `None`.
- `load_result` reconstructs QP or DMRG `Eigenstates` from the JSON and its
  companion wavefunction file. It also reads the earlier QP result format.

The modules `qdjj_solver.qp_solver` and `qdjj_solver.dmrg_solver` also provide
their own `solve` and `SolverOptions`. When calling either module's `solve`
directly, use its `SolverOptions`; the top-level `qdjj_solver.solve` additionally
accepts a dictionary. Starting states are described [below](#starting-states-and-continuation).

### Defaults

| Setting | Python | `qdjj-solver` command |
|---|---|---|
| Numerical method | QP | QP |
| QP cutoff | `None`: all bath configurations allowed by the sector | `2` |
| Reference-model compression | `reference_model(..., compress=True)` | `false` |
| Sector | Odd parity, `twice_sz=1`, unrestricted eta/number | Same |
| Requested eigenstates | `1` | `1` |

Compression applies only to eligible symmetric reference models. Use
`compress=False` explicitly in Python when the full finite-bath spectrum is
intended. DMRG always requires `cutoff=None`, and its default
`require_convergence=False` returns diagnostics even for coarse, unconverged calculations.
Set `require_convergence=True` for calculations that must meet all convergence
checks for the finite Hamiltonian. QP results must always pass the projected-residual and orthonormality checks.

### QP solver options

Use `qdjj_solver.qp_solver.SolverOptions` or a dictionary of its fields. All fields, defaults and
units are listed in the [option reference](reference/options.md#qp-solver).

Requesting nearly the whole spectrum can force dense
diagonalization even with an iterative method.

### DMRG solver options

Use `qdjj_solver.dmrg_solver.SolverOptions`; this requires installation with `[dmrg]`. All fields,
defaults and units are listed in the [option reference](reference/options.md#dmrg-solver).

Require `min_sweeps <= max_sweeps` and, with the mixer enabled,
`mixer_sweeps < min_sweeps`. `require_convergence=True` requires residual
evaluation. Disabling residuals leaves `finite_problem_converged=False`.
Residual evaluation has its own temporary memory cost, described in
[numerics](numerics.md#residuals-and-convergence-checks).

The local Lanczos threshold `lanczos_probability_tolerance` (TeNPy's `P_tol`, default
`1e-14`) is independent of `energy_tolerance` and of the residual tolerance
`residual_tolerance`. Check the full-state residual explicitly to assess the combined accuracy. The
solver passes the local threshold to TeNPy for each calculation, preserving TeNPy's defaults. Saved results
record the value used in `metadata["options"]`.

### Results and transitions

The fields of both result types are listed in the [result reference](reference/results.md):
`energies`, `residuals`, `observables`, `hamiltonian`, `metadata` and `timings` for both
methods; `vectors`, `basis` and `qp_weights` for QP; `states` and `prepared` for DMRG. State indices
are zero-based.

Named observable definitions and units are in the [observable reference](reference/observables.md).
For QP, all named observables are measured by default; the backend's `solve` and
`QPWorkspace.solve` accept `observables` to select a subset. DMRG's `observables`
option can also select a subset. Residuals measure how well the state satisfies the eigenvalue equation
for the finite Hamiltonian; QP residuals refer to the chosen projected space.

The snippets below continue this chapter's `h`, `sector`, `qp`, and (with the
optional DMRG installation) `mps` calculation. **Functions on DMRG results** return
overlaps or transition matrices with bra eigenstates as rows and ket eigenstates as columns:

```python
gram = mps.overlaps()
matrix = mps.matrix_elements(h.observables["impurity_spin_z"])
```

`bra_result.matrix_elements(operator, ket_result)[i, j]` evaluates
$\langle\mathrm{bra}_i\vert\mathrm{operator}\vert\mathrm{ket}_j\rangle$.
DMRG transition results require identical canonical coordinates, mode ordering,
grouping, and types of conserved quantum numbers. Their values may differ. For example,
both results must use parity-only tensors or both must use parity-plus-spin tensors.

**QP results** offer the same functions. The operator may be a `FermionOperator`
in canonical coordinates or a probe name (`h.probe(name)`; see the
[probe list](reference/observables.md#probes)); any other name selects a model
observable. The ket result may belong to another sector or QP cutoff.

```python
odd = solve(h, cutoff=2, sector=Sector(1, 1))
even = solve(h, cutoff=2, sector=Sector(0, 0), options={"eigenpairs": 2})
addition = odd.matrix_elements("d_up_dagger", even)   # <D|d_up^dagger|S>, <D|d_up^dagger|S*>
pair = even.matrix_elements("current")                 # includes <S|dH/dphi|S*>
```

Requesting matrix elements between sectors that no term of the operator connects raises an error. For
single vectors, `ProjectedOperator` gives direct access:

```python
from qdjj_solver.qp_solver import ProjectedOperator

spin = ProjectedOperator(h.observables["impurity_spin_z"], qp.basis)
mean_spin = spin.expectation(qp.vectors[:, 0])
```

`ProjectedOperator(operator, ket_result.basis).transition(bra, bra_basis, ket)`
acts between compatible Fock bases, allowing different sectors and cutoffs.
The caller is responsible for supplying operators and vectors in the same canonical
coordinates, in addition to satisfying the dimension and spin-label checks.
Use `h.physical_operator(electron_operator)` to transform physical electron
operators through an uncompressed model's electron-to-canonical transformation.

### Local orbital density matrices

`orbital_density_matrix` works with both QP and DMRG eigenstates. Orbital and
state indices are zero-based. Its basis is **empty, up, down, up-down**, with
the doubly occupied state defined by creating up before down.

```python
from qdjj_solver import orbital_density_matrix, orbital_density_properties

rho = orbital_density_matrix(qp, orbital=0, state=0)  # shape (4, 4)
properties = orbital_density_properties(rho)
print(properties["empty"], properties["single"], properties["double"])
print(properties["pair_amplitude"], properties["plus"], properties["minus"])
all_rho = orbital_density_matrix(qp, orbital=0, state=None)  # shape (k, 4, 4)
# Equal mixture of a complete two-state multiplet, for a result with k=2:
# mixed = orbital_density_matrix(qp, 0, state=None, weights=[1, 1])
```

`weights` must be real, finite, nonnegative, one per retained state, and have a
positive sum. They are normalized to sum one and describe an **incoherent**
mixture. A complete, equally weighted degenerate eigenspace has a
basis-independent RDM; an individual degenerate eigenvector generally does not.

The default `physical=True` transforms physical impurity-electron operators
through the saved coordinate map. The impurity-only eta map also supports
verified paired-mode compression. For a manually constructed canonical-mode
Hamiltonian without an electron map, pass `physical=False` explicitly.
This API concerns impurity orbitals, not arbitrary reservoir subsystems.

`orbital_density_operators(orbital=0, hamiltonian=h)` supplies the canonical
operators for independent measurement or matrix caching. It returns a mapping
`(row,column) -> operator`, where the expectation equals `rho[row,column]`.
Without `hamiltonian`, it returns physical-electron expressions. Composite
projectors are built before QP truncation, preserving vacuum contractions.
See [density-matrix conventions](models.md#local-orbital-density-matrix) for
the signs, coherences, and paired-state probabilities.

### Floquet calculations

Periodic drives of model parameters are described by `PeriodicHamiltonian`,
which samples a model builder over one period:

```python
from qdjj_solver import ParameterDrive, PeriodicHamiltonian, reference_model
from qdjj_solver.qp_solver import floquet_propagator, rwa_model, solve_floquet

periodic = PeriodicHamiltonian.from_parameters(
    reference_model, dict(bath=bath, phi=.8, symmetry=False, compress=False),
    ParameterDrive.cosine("phi", .05), omega=.57)
floquet = solve_floquet(periodic, cutoff=2, sector=sector, options=dict(floquet_cutoff=4, states=2))
```

`solve_floquet`, `floquet_propagator` and `rwa_model` return `FloquetStates`
with quasienergies, harmonic weights, static-state weights and Floquet modes.
`periodic.van_vleck(order)` returns a static high-frequency `Hamiltonian` for
either backend. Options are collected in `FloquetOptions`
([option reference](reference/options.md#floquet-calculations)). Persistence through
`save_result` currently supports static QP and DMRG eigenstates. Conventions, methods, result fields and
convergence checks are described in [Floquet calculations](floquet.md).

### Starting states and continuation

QP iterative solvers accept a finite, nonzero vector of shape `(D,)` in the
target basis ordering. For example, refine the same sector and cutoff with:

```python
continued_qp = solve(h, cutoff=2, sector=sector,
                     options={"method": "sparse", "tolerance": 1e-12},
                     initial=qp.vectors[:, 0])
```

The vector is used by sparse/matrix-free iterations. The block eigensolver also
accepts a starting array of shape `(D, n)` and supplements it with random columns
up to the requested number of states. Dense and Schur methods compute their
states directly. QP continuation across cutoffs requires an explicit embedding
into the new Fock basis. A [QP workspace](qp_workflow.md#reuse-preparation-and-check-spectrum-coverage)
reuses compatible preparation and starting states across repeated calculations.

DMRG takes a prior DMRG `Eigenstates` object:

```python
continued_mps = solve(h, sector=sector, backend="dmrg", initial=mps,
                      options={"eigenpairs": 2, "chi_max": 64,
                               "require_convergence": True})
```

The canonical basis, MPS mode ordering and grouping, types of conserved quantum
numbers, and sector must match. Parameter changes are allowed when the
electron-to-canonical transformation remains identical; changing bath nodes
generally changes that transformation. The first trial for each previously
calculated eigenstate uses its previous MPS; other trials remain independent.

## Input files and terminal commands

The `qdjj-solver` command reads a JSON input file, with settings enclosed in
`{...}`. Version 1 accepts the following entries:

| Field | JSON type | When omitted |
|---|---|---|
| `format_version` | Integer `1` | `1` |
| `backend` | String `"qp"` or `"dmrg"` | `"qp"` |
| `cutoff` | Integer from `0` through `h.bath_modes`, or `null` for full space | `2` for QP; `null` for DMRG |
| `model` | Object | [Reference-model defaults](models.md#reference_model), with `compress=false` in `qdjj-solver` |
| `bath` | Object | `{"pairs": 8}` for a reference model; other model kinds specify their own modes/reservoirs |
| `sector` | Object | `{"parity": 1, "twice_sz": 1}`; unrestricted eta/number |
| `solver` | Object | Chosen method's `SolverOptions` defaults |

Unknown setting names are reported with their location, such as
`model.impurity.terms[0]`, to help identify errors in the input. Use JSON
`true`/`false` for yes/no settings and integers for the cutoff and sector labels.
Where a bath or model specification accepts `metadata`, that dictionary can
contain your own notes. Place such notes inside the supported bath or model entry.
See the [examples](../examples/README.md).
The command-line option `--backend` overrides the `backend` of the input file. DMRG requires
an omitted or null `cutoff`. Cutoff and solver options must be valid for the selected method,
also when `--backend` overrides it.
The `qdjj-solver` command defaults to QP with cutoff two; Python `solve` defaults instead
to full finite space (`cutoff=None`). DMRG `--resume` names a result JSON with
a companion saved MPS file.

`--resume` initializes a new optimization from the previously saved states;
settings such as `chi_max` and tolerances come from the new configuration. The
saved file supplies completed MPS states as the initial wavefunctions for this
new optimization. The canonical basis, MPS mode ordering and grouping,
sector, and types of conserved quantum numbers must be compatible.
See [`examples/dmrg_refined.json`](../examples/dmrg_refined.json) for a refinement
of [`examples/dmrg.json`](../examples/dmrg.json).

The command uses an uncompressed reference-model bath by default for both
methods. Switching only the method therefore preserves the finite Hamiltonian;
any dark-mode restriction must be requested explicitly with `compress=true`.

Model kinds:

- `reference`: parameters for [`reference_model`](models.md#reference_model), plus a top-level bath.
- `multi-orbital`: `impurity` with `orbitals` and polynomial `terms`, `reservoirs`,
  `tunneling`, and optional `direct` entries with `leads` and `matrix`. An optional
  `direct_derivatives` entry associates observable names with lists in the same contact
  format. Duplicate direct contacts are rejected. Other [`make_model`](models.md#make_model)
  options, including `phases`, `phase_velocities`, `bath_reference`, `eta_basis`,
  `compress_pairs`, and `coefficient_tolerance`, also belong inside `model`.
- `hamiltonian`: the fields of `Hamiltonian.record()`, plus `kind`.

Bath kinds are `cosh-grid` (the default; the arguments of `cosh_grid`), `surrogate`
(the arguments of `fit_surrogate`), `chain-expansion` (`levels`, `delta`, `bandwidth` and
the optional boolean `wide_band`), and `discrete` (`xi`, `weights`, `delta`, `bandwidth` and
optional `metadata`). Complex arrays use `{"real": [...], "imag": [...]}`. Impurity
observables may also be given as lists of fermion-operator terms in a multiorbital input
file.

For example, this complete QP input describes two interacting orbitals and two
reservoirs. Its literal impurity terms specify the energy zero; include the
constant explicitly to reproduce the convention of `Impurity.anderson`. Operator indices in `terms` are signed
and one-based, while reservoir indices in `leads` are zero-based:

```json
{
  "format_version": 1,
  "backend": "qp",
  "model": {
    "kind": "multi-orbital",
    "impurity": {
      "orbitals": 2,
      "terms": [
        {"operators": [1, -1], "real": -0.4, "imag": 0.0},
        {"operators": [2, -2], "real": -0.4, "imag": 0.0},
        {"operators": [3, -3], "real": -0.3, "imag": 0.0},
        {"operators": [4, -4], "real": -0.3, "imag": 0.0},
        {"operators": [1, 2, -2, -1], "real": 0.8, "imag": 0.0},
        {"operators": [3, 4, -4, -3], "real": 0.6, "imag": 0.0}
      ]
    },
    "reservoirs": [
      {"kind": "cosh-grid", "pairs": 1, "delta": 1.0, "bandwidth": 5.0},
      {"kind": "cosh-grid", "pairs": 1, "delta": 1.0, "bandwidth": 5.0}
    ],
    "tunneling": [
      [[0.2, 0.0, 0.0, 0.0], [0.0, 0.2, 0.0, 0.0]],
      [[0.0, 0.0, 0.15, 0.0], [0.0, 0.0, 0.0, 0.15]]
    ],
    "phases": [-0.4, 0.4],
    "phase_velocities": [-0.5, 0.5],
    "direct": [{"leads": [0, 1], "matrix": [[0.02, 0.0], [0.0, 0.02]]}],
    "bath_reference": "isolated",
    "eta_basis": false,
    "compress_pairs": false,
    "coefficient_tolerance": 1e-14
  },
  "cutoff": 2,
  "sector": {"parity": 0, "twice_sz": 0},
  "solver": {"eigenpairs": 1}
}
```

The [bath representation guide](bath_representations.md#switching-in-python-and-configuration-files)
provides examples of changing bath representations and reusing saved fits, all generator
parameters, and accuracy/time benchmarks.
The [ChE guide](chain_expansion.md) explains its low-frequency Padé construction
and the density-normalization meaning of `bandwidth` when `wide_band=true`.

## Saved results and reproducibility

For a complete scan, scalar-only analysis and saved-wavefunction workflow, see
[Parameter scans, saved results and reproducible analysis](workflows.md).
The [record reference](reference/records.md) specifies nested fields and companion files.

The JSON format `qdjj-eigenstates`, version 1, contains the numerical method,
energies, residuals, observables, elapsed times, solver settings, symmetry
sectors, Hilbert-space restrictions, and convergence diagnostics. It also saves
the complete Hamiltonian, software versions, and SHA-256 checksums identifying
the Hamiltonian and the solver files used. Wavefunctions are optional and saved
in a separate file. Non-finite JSON numbers are rejected; undefined initial
sweep differences are stored as `null`.
The Python result's `metadata` is stored under `calculation` in JSON, and
`timings` under `timings_seconds`.

`qdjj_solver.common.io.source_provenance()` calculates checksums of the Python modules and
C++ library actually used by the calculation. It stores their package-relative
names, so the result can be compared between computers with different
installation directories.

```python
from qdjj_solver import save_result, load_result
save_result(qp, "results/state.json", save_states=True)
restored = load_result("results/state.json")
```

`qdjj_solver.save_result` uses this format for both methods. `load_result`
identifies the method from the saved file and requires the companion
wavefunctions. Results saved without wavefunctions can still be read with
Python's standard `json` module. `restored.backend` is `"qp"` or `"dmrg"`.
The base installation suffices for loading QP states.
To move a calculation, copy its JSON and the companion file named in
`record["state_artifact"]["file"]` together into the destination directory.

The DMRG wavefunction file `state.<unique-id>.mps.h5` is identified by
`kind="tenpy-mps-hdf5"`. Its checksum and recorded basis identify the matching
JSON result and prevent combining files from different saves. Reading the MPS
depends on TeNPy's support for its saved-state format; the JSON remains readable
independently. Use TeNPy MPS files from trusted calculations.

QP wavefunction files have `kind="fock-amplitudes-npz"` and contain energies,
eigenvector columns, Fock occupations, and QP counts. They are read with NumPy's
`allow_pickle=False`. New NPZ files identify the matching JSON by its
`result_sha256` checksum. Loading verifies both file checksums, the energies,
basis ordering, dimensions, and orthonormality. Earlier four-array NPZ files
remain readable, with checksum, energy, and basis checks.

Always obtain the companion filename from the JSON record. Each save uses a
unique wavefunction filename and replaces the JSON only after that file has
been written successfully. A failed save therefore preserves the previous
result. Older wavefunction files remain available for copies of the earlier
JSON. They can be removed once they are no longer needed by a retained result
or a running analysis. An interrupted save may also leave an unused wavefunction
file.

Version 1 allows additional information to be saved. Loading rejects unsupported
JSON or wavefunction-format versions. Reading MPS files also depends on the
installed TeNPy version.

## Troubleshooting and exit status

Exit codes, common error messages and their remedies are collected in
[troubleshooting](troubleshooting.md).
