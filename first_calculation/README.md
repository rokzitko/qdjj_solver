# First calculation: inputs and reference results

This directory reproduces the [main README's first calculation](../README.md#first-calculation),
including its optional DMRG and exact-diagonalization comparison.

## Run

Run these commands from the repository root after installing the solver:

```sh
python -m pip install .
python -B first_calculation/run.py
```

The default run performs the initial QP calculation, checks the README's energy,
residual, and normalization assertions, and saves `results/first_calculation/qp.json`.
To include the two-state DMRG and exact calculations:

```sh
python -m pip install '.[dmrg]'
python -B first_calculation/run.py --all
```

The full run checks that all three calculations use the same Hamiltonian and
sector, checks DMRG against exact energies, and saves `qp.json`, `dmrg.json`,
`exact.json`, and `comparison.json` under `results/first_calculation/`.
Use `--output DIRECTORY` to choose another destination. Input paths are resolved
relative to `run.py`, so the script can also be invoked from another directory.

## Inputs and physical conventions

| Input | Calculation | Requested states | Reference output |
|---|---|---|---|
| [`input/qp.json`](input/qp.json) | QP cutoff two | 1 | [`reference/qp.json`](reference/qp.json) |
| [`input/dmrg.json`](input/dmrg.json) | DMRG, `chi_max=32`, convergence required | 2 | [`reference/dmrg.json`](reference/dmrg.json) |
| [`input/exact.json`](input/exact.json) | QP backend, `cutoff=null` (full finite space) | 2 | [`reference/exact.json`](reference/exact.json) |

All inputs specify the same uncompressed, symmetric two-lead reference junction:

- `cosh-grid`, `pairs=1`, gap $\Delta=1$, normal-state **half-bandwidth** $D=5$;
- interaction $U=1.6$, total hybridization $\Gamma=0.3$ ($0.15$ per lead);
- phase difference $\phi=0.8$ radians, with lead phases $[-\phi/2,+\phi/2]$;
- default zero detuning and zero field;
- odd parity and `twice_sz=1`, namely $S_z=1/2$;
- `compress=false`, seed `1729`, and one numerical thread.

Energies use the gap as the unit and the solver's common
[energy reference](../docs/models.md#units-and-energy-reference).
The QP input requests one state to match the first Python snippet; the separate
[`examples/qp.json`](../examples/qp.json) command-line example requests two.

Each input also works directly with the unified CLI:

```sh
qdjj-solver solve first_calculation/input/qp.json --output results/first_calculation/qp.json
qdjj-solver solve first_calculation/input/dmrg.json --output results/first_calculation/dmrg.json
qdjj-solver solve first_calculation/input/exact.json --output results/first_calculation/exact.json
```

The runner additionally performs the cross-method checks and writes the
comparison file.

## Reference outputs and checks

The three result files use the standard `qdjj-eigenstates` JSON format. They
contain energies, residuals, observables, numerical settings, the complete
Hamiltonian and bath records, software versions, and implementation checksums.
The QP and exact results also include QP-number probabilities.

The archived energies, in gap units, are:

| Calculation | State 0 | State 1 |
|---|---|---|
| QP cutoff two | `-0.295059596719` | — |
| DMRG | `-0.296312491144` | `1.727646371364` |
| Exact finite bath | `-0.296312491144` | `1.727646371364` |

The initial QP result has QP weights approximately
`[0.89270561, 0.10104019, 0.00625420]` for zero, one, and two bath quasiparticles.
Its lowest energy is about `0.001252894425` above the exact value. The largest
DMRG/exact energy difference in the archived run is below `9e-16`.
The runner checks:

- absolute energy error against the README value below `1e-10`;
- projected QP residual below `1e-9`;
- QP-weight normalization error below `1e-12`;
- DMRG finite-problem convergence and DMRG/exact energy agreement with
  `rtol=0`, `atol=1e-9` when `--all` is selected.

[`reference/comparison.json`](reference/comparison.json) records the DMRG-minus-exact
energy differences, the cutoff-two shift of the lowest energy relative to exact,
and the DMRG impurity-spin matrix printed by the README. Matrix rows and columns
are bra and ket state indices, respectively, ordered by energy starting at zero;
spin is in units of $\hbar$. Complex arrays use `real` and `imag` fields.
Eigenstate phases can change off-diagonal matrix elements between runs, and
degenerate eigenstates may rotate within their common subspace.

Use numerical tolerances when comparing a fresh run with the references.
Residuals are bounded accuracy checks; timings, software versions, implementation
checksums, and last-place floating-point digits can differ between environments.
The selected sector's lowest state does not establish the global ground state;
compare the relevant even and odd sectors for that. QP cutoff/MPS bond convergence
and bath convergence are separate from the eigensolver checks above.

### Saving and loading wavefunctions

Add `--save-states` to the runner or to any of the CLI commands to save companion
NPZ (QP/exact) or HDF5 (DMRG) files:

```sh
python -B first_calculation/run.py --all --save-states
```

Then, for example:

```python
from qdjj_solver import load_result

qp = load_result("results/first_calculation/qp.json")
print(qp.energies, qp.qp_weights)
```

The checked-in reference files contain scalar results and can be read with
Python's `json` module. `load_result` requires the companion wavefunctions from
a run using `--save-states`; keep each JSON and its named companion together.

## Regenerate the references

From the repository root, with the DMRG extra installed:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B first_calculation/run.py --all --output first_calculation/reference
```

This runs the numerical checks before saving the reference files. Review the
numerical changes and recorded provenance when updating them. Routine runs use
the ignored `results/first_calculation/` destination.
