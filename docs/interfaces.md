# Using the solver and saving results

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
Set `require_convergence=True` for calculations that must meet all finite-problem
checks. QP results must always pass the projected-residual and orthonormality checks.

### QP solver options

Use `qdjj_solver.qp_solver.SolverOptions`. Energy units below are those of the
input Hamiltonian; iteration counts and probabilities are dimensionless.

| Option | Default | Meaning and units |
|---|---|---|
| `method` | `"auto"` | `"dense"`, `"sparse"`, `"matrix-free"`, or `"schur"`; see [numerics](numerics.md#qp-diagonalization) |
| `eigenpairs` | `1` | Number of lowest eigenstates in the selected sector |
| `tolerance` | `1e-11` | Relative ARPACK eigenvalue tolerance; Schur root-search tolerance in energy units, tightened internally near its pole; unused by dense diagonalization |
| `residual_tolerance` | `2e-9` | Maximum projected residual norm, in energy units |
| `maxiter` | `10000` | Maximum ARPACK iterations |
| `ncv` | `40` | Requested Krylov-space size; increased as needed for the number of eigenstates and limited by the sector dimension |
| `max_dimension` | `5_000_000` | Maximum number of Fock basis states, checked before enumeration |
| `max_nnz` | `100_000_000` | Entry limit for sparse matrix construction, used by dense, sparse, and Schur methods; the remaining `max_memory_gib` budget can lower this limit |
| `max_memory_gib` | `24.0` | Limit on estimated basis/matrix storage in GiB; actual total memory use can be larger |
| `threads` | `1` | Maximum number of computational threads used by the numerical libraries during diagonalization/refinement |
| `seed` | `1729` | Random-number seed for an automatically generated iterative starting vector |

`auto` selects dense for dimension at most 128, sparse up to 400,000, and
matrix-free above that. Requesting nearly the whole spectrum can force dense
diagonalization even with an iterative method. Every returned QP result passes
the projected residual and orthonormality checks.

### DMRG solver options

Use `qdjj_solver.dmrg_solver.SolverOptions`; this requires installation with `[dmrg]`.

| Option | Default | Meaning and units |
|---|---|---|
| `eigenpairs` | `1` | Number of eigenstates sought sequentially |
| `chi_max` | `128` | Maximum MPS bond dimension |
| `max_sweeps` | `60` | Maximum sweeps per initial-state trial for each eigenstate |
| `min_sweeps` | `10` | Minimum sweeps per initial-state trial for each eigenstate |
| `energy_tolerance` | `1e-11` | Absolute sweep-energy stability scale, in energy units; also controls local Lanczos `E_tol` |
| `entropy_tolerance` | `1e-8` | Dimensionless entropy-change threshold in TeNPy's stopping test |
| `residual_tolerance` | `1e-7` | Maximum physical residual norm for finite-problem convergence, in energy units |
| `orthogonality_tolerance` | `1e-8` | Dimensionless tolerance on deviations of the overlap matrix from the identity and overlap with earlier eigenstates |
| `svd_min` | `1e-14` | Minimum retained singular value in optimization; dimensionless, may be zero |
| `lanczos_maxiter` | `100` | Maximum local Lanczos iterations |
| `lanczos_probability_tolerance` | `1e-14` | Positive, finite, dimensionless local Lanczos probability-error threshold; passed to TeNPy as `P_tol` |
| `seed_trials` | `4` | Number of initial-state trials per eigenstate, including independently randomized states |
| `seed` | `1729` | Random-number seed for local transformations preserving the chosen quantum numbers |
| `seed_randomization_steps` | `4` | Number of layers of local random unitary transformations for additional initial states; zero disables randomization |
| `excitation_operator` | `None` | Name of a model observable that preserves the chosen quantum numbers and generates an initial excited state |
| `mixer` | `True` | Enable TeNPy's optimization mixer |
| `mixer_amplitude` | `1e-5` | Initial mixer strength |
| `mixer_sweeps` | `6` | Sweep at which the mixer is disabled |
| `group_size` | `2` | One or two canonical fermionic modes per MPS site |
| `mode_order` | `None` | Permutation of all zero-based canonical mode indices; `None` uses their order in the model |
| `observables` | `None` | Measure all named observables; a tuple selects names, and `()` selects none |
| `calculate_residuals` | `True` | Evaluate the uncompressed physical residual after optimization |
| `require_convergence` | `False` | Stop with an error if the finite-problem convergence checks fail |
| `threads` | `1` | Numerical-library thread limit during optimization and measurements |

Require `min_sweeps <= max_sweeps` and, with the mixer enabled,
`mixer_sweeps < min_sweeps`. `require_convergence=True` requires residual
evaluation. Disabling residuals leaves `finite_problem_converged=False`.
Residual evaluation has its own temporary memory cost, described in
[numerics](numerics.md#residuals-and-convergence-checks).

The probability tolerance is independent of `energy_tolerance` and the global
`residual_tolerance`; there is no automatic mapping that guarantees a global residual.
The `1e-14` default is unchanged; no upstream TeNPy defaults are modified.
Saved `metadata["options"]` records the actual `lanczos_probability_tolerance`.

### Results and transitions

Let $k$ be the number of calculated eigenstates, $D$ the QP basis dimension, and $q$ the
effective QP cutoff (`h.bath_modes` when `cutoff=None`). State indices are zero-based.

| Field | Method | Array dimensions or content; interpretation |
|---|---|---|
| `backend` | Both | `"qp"` or `"dmrg"` |
| `energies` | Both | `(k,)`; eigenenergies in Hamiltonian units |
| `residuals` | Both | `(k,)`; projected norms for QP, physical norms for DMRG; `None` if DMRG residuals were disabled |
| `observables` | Both | Dictionary of `(k,)` expectation arrays; values may be complex for non-Hermitian operators |
| `hamiltonian` | Both | The shared `Hamiltonian` used in the calculation |
| `metadata` | Both | Numerical settings, symmetry sector, Hilbert-space restrictions, and convergence diagnostics |
| `timings` | Both | Dictionary of elapsed times in seconds |
| `vectors` | QP | `(D, k)`; eigenvector `i` is column `vectors[:, i]` |
| `basis` | QP | `FockBasis`; defines vector-row ordering and canonical occupations |
| `qp_weights` | QP | `(k, q+1)`; `qp_weights[i, n]` is the probability of `n` bath QPs in state `i` |
| `states` | DMRG | List of $k$ TeNPy MPS objects |
| `prepared` | DMRG | Canonical basis, conserved quantum numbers, and MPS mode ordering and grouping |

Named observable definitions and units are in [models](models.md#built-in-observables).
For QP, all named observables are measured. DMRG's `observables` option can select
a subset. Residuals measure how well the state satisfies the eigenvalue equation
for the finite Hamiltonian; QP residuals refer to the chosen projected space.

The snippets below continue the README's `h`, `sector`, `qp`, and (with the
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
parity-only tensors cannot directly contract with parity-plus-spin tensors.

For **QP results**, use `ProjectedOperator`:

```python
from qdjj_solver.qp_solver import ProjectedOperator

spin = ProjectedOperator(h.observables["impurity_spin_z"], qp.basis)
mean_spin = spin.expectation(qp.vectors[:, 0])
```

`ProjectedOperator(operator, ket_result.basis).transition(bra, bra_basis, ket)`
acts between compatible Fock bases, allowing different sectors and cutoffs.
The caller supplies operators and vectors in the same canonical coordinates;
dimension and spin-label checks do not identify an arbitrary change of physical
basis. Use `h.physical_operator(electron_operator)` to transform physical electron
operators through an uncompressed model's electron-to-canonical transformation.

### Starting states and continuation

QP iterative solvers accept a finite, nonzero vector of shape `(D,)` in the
target basis ordering. For example, refine the same sector and cutoff with:

```python
continued_qp = solve(h, cutoff=2, sector=sector,
                     options={"method": "sparse", "tolerance": 1e-12},
                     initial=qp.vectors[:, 0])
```

The vector is used by sparse/matrix-free iterations. Dense and Schur methods
compute their states directly. QP continuation across cutoffs requires an
explicit embedding into the new Fock basis.

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
contain your own notes; `metadata` is not a top-level calculation setting.
See the [examples](../examples/README.md).
`--backend` selects the method in preference to the input file. DMRG requires omitted or null `cutoff`.
Cutoff and solver options must be valid for the selected method, including when
using this override.
The `qdjj-solver` command defaults to QP with cutoff two; Python `solve` defaults instead
to full finite space (`cutoff=None`). DMRG `--resume` names a result JSON with
a companion saved MPS file.

`--resume` initializes a new optimization from the previously saved states;
settings such as `chi_max` and tolerances come from the new configuration. The
saved file contains completed MPS states; it cannot restart an interrupted sweep
at its stopping point. The canonical basis, MPS mode ordering and grouping,
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

Bath kinds are `cosh-grid`, `surrogate` (fit settings), `chain-expansion`
(`levels`, `delta`, `bandwidth`, optional boolean `wide_band`), and `discrete`
(explicit xi, weights, gap, bandwidth, optional metadata). Complex arrays use
`{"real": [...], "imag": [...]}`. Impurity observables may also be given as
lists of fermion-operator terms in a multiorbital input file.

For example, this complete QP input describes two interacting orbitals and two
reservoirs. Its literal impurity terms specify the energy zero; they do not add
the constant used by `Impurity.anderson`. Operator indices in `terms` are signed
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

The JSON format `qdjj-eigenstates`, version 1, contains the numerical method,
energies, residuals, observables, elapsed times, solver settings, symmetry
sectors, Hilbert-space restrictions, and convergence diagnostics. It also saves
the complete Hamiltonian, software versions, and SHA-256 checksums identifying
the Hamiltonian and the solver files used. Wavefunctions are optional and saved
in a separate file. Non-finite JSON numbers are rejected; undefined initial
sweep differences are stored as `null`.
The Python result's `metadata` is stored under `calculation` in JSON, and
`timings` under `timings_seconds`.

`common.io.source_provenance()` calculates checksums of the Python modules and
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
Loading QP states does not require TeNPy.
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

The `qdjj-solver solve` command returns exit status **0** after saving a result
and **2** for recognized errors in input, required libraries, numerical calculations,
or file access. Invalid command arguments also return 2. With DMRG's default
`require_convergence=False`, exit 0 can
accompany `Finite-problem checks passed: False`; use `require_convergence=true`
in JSON when a calculation must meet those checks. Bath convergence is assessed
separately even when finite-problem checks pass.

| Message or symptom | What to check |
|---|---|
| DMRG/TeNPy or HDF5 dependency missing | From the source directory, run `python -m pip install '.[dmrg]'` with the interpreter used to run the calculation. |
| `configuration` / `model` / `solver` must be an object | Use `{...}` for these objects; use JSON booleans and integer sector values. Check the reported nested field path for typos. |
| Hamiltonian does not conserve `S_z` | For spin-mixing terms choose `Sector(parity=..., twice_sz=None)`, or JSON `"twice_sz": null`. |
| Fermion parity and `twice_sz` are incompatible | Even parity requires an even integer `twice_sz`; odd parity requires an odd integer. |
| Empty symmetry sector / more eigenpairs than states | Check quantum numbers, bath cutoff, included modes, and number of requested eigenstates. |
| Invalid QP cutoff | Use `0 <= cutoff <= h.bath_modes` after model compression, or `None`/JSON `null` for full space. Impurity-only QP CLI inputs need explicit `cutoff: 0` or `null`, not the default `2`. |
| Sector dimension or memory estimate exceeds a limit | Estimate the basis size before increasing `max_dimension` or `max_memory_gib`. Consider a converged QP cutoff or a DMRG calculation. |
| Sparse matrix exceeds `max_nnz` | Dense, sparse, and Schur methods build this matrix. Check both `max_nnz` and the memory-derived cap from `max_memory_gib`; increasing only `max_nnz` may not help. Consider `method="matrix-free"`; its basis and Krylov-space limits still apply. |
| DMRG requires `cutoff=null` | Omit the JSON cutoff or set it to `null`. MPS bond dimension is controlled by `chi_max`. |
| Finite-problem convergence failed | Inspect residuals, orthonormality, energy ordering, and sweep histories. Increase bond dimension, sweeps, or the number of initial-state trials as appropriate. |
| Initial states have incompatible coordinates/layout | Reusing MPS states requires the same canonical basis, mode ordering, and site grouping. Start a new calculation if these change. |
| No supported state artifact / checkpoint missing | Save with `save_states=True` or `--save-states`, then retain both the JSON and its named wavefunction file. |
| Checksum or scalar-record mismatch | Recover the matching JSON and wavefunction files from the same save. Editing the JSON also invalidates the checksum linking it to a new-format wavefunction file. |

Consult [numerics](numerics.md) for convergence choices. To report a reproducible
issue, include a small configuration, the error, the software version, and the
result's saved software-version information when available. Report problems on
the [project's issue page](https://github.com/rokzitko/qdjj_solver/issues).
