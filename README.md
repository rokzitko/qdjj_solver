# qdjj_solver

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23104978.svg)](https://doi.org/10.5281/zenodo.23104978)

**Quasiparticle (QP) expansion and density-matrix renormalization group (DMRG)
solvers for superconducting quantum impurities.**

Authors: **Teodor Iličin and Rok Žitko**. License:
[BSD-3-Clause](https://github.com/rokzitko/qdjj_solver/blob/main/LICENSE).

`qdjj_solver` calculates low-lying many-body energies, impurity charge and spin,
Josephson currents, and transition matrix elements for interacting impurities
coupled to superconducting reservoirs. Models can contain one or more impurity
orbitals and several superconducting leads.

Both numerical methods solve the same finite-bath Hamiltonian:

| Method | Approximation and convergence parameter |
|---|---|
| Quasiparticle expansion | Retains all impurity configurations and limits the number of bath QPs. Increase the QP cutoff to check convergence; removing it gives exact diagonalization (ED) of the finite-bath model in the selected symmetry sector. |
| DMRG | Represents the state as a matrix product state (MPS). Increase the MPS bond dimension and check the optimization to establish convergence for the finite bath. |

Bath discretization is a separate approximation. Small finite models allow
direct comparisons between the two methods and unrestricted ED.

The [QP method references](https://github.com/rokzitko/qdjj_solver/blob/main/docs/references.md#qp-variational-method)
describe our Andreev spin-qubit (ASQ) variational theory and its earlier single-impurity and
quantum-dot Josephson-junction formulations.

The Python modules `qdjj_solver.common`, `qdjj_solver.qp_solver`, and
`qdjj_solver.dmrg_solver` provide the common model definitions and the two solvers.

## Install

Currently, install from [source](https://github.com/rokzitko/qdjj_solver).
Building the native extension requires **Python 3.11+**, a **C++17-capable
compiler**, and **Python development headers** for the interpreter being used.
Install the platform toolchain before running pip:

- **Linux:** GCC or Clang; on Debian/Ubuntu, for example, install
  `build-essential` and `python3-dev` (or the development package matching your
  Python version).
- **macOS:** Xcode Command Line Tools (`xcode-select --install`).
- **Windows:** Visual Studio C++ Build Tools with the C++ workload and Windows SDK.

Pip obtains CMake and pybind11 as needed, but does **not** install a C++ compiler
or Python development headers.

```sh
git clone https://github.com/rokzitko/qdjj_solver.git
cd qdjj_solver
```

Choose **one** of these installation alternatives. The **`[dmrg]` extra is
required for every DMRG example below**, including continuation and
`first_calculation/run.py --all`; it installs TeNPy and HDF5 support for saved
MPS states. The base installation supports models and the QP solver without TeNPy.

```sh
python -m pip install .                 # models and QP solver
python -m pip install '.[dmrg]'         # add DMRG and saving of MPS states (recommended)
python -m pip install -e '.[dmrg,dev]'  # development installation
```

## First calculation

```python
from qdjj_solver import Sector, cosh_grid, reference_model, solve

bath = cosh_grid(pairs=1, delta=1.0, bandwidth=5.0)
h = reference_model(bath, u=1.6, gamma=0.3, phi=0.8, compress=False)
sector = Sector(parity=1, twice_sz=1)  # odd parity, S_z = 1/2
qp = solve(h, cutoff=2, sector=sector, backend="qp")
print(qp.energies, qp.qp_weights)

assert abs(qp.energies[0] - (-0.2950595967187083)) < 1e-10
assert max(qp.residuals) < 1e-9
assert abs(qp.qp_weights[0].sum() - 1) < 1e-12
```

From the main source directory, the same QP-only example runs with
`python -B examples/quickstart.py`.
It computes the lowest state in the selected sector of a small finite model.
Here energies are in units of the gap, `gamma=0.3` is the total hybridization
of both leads, and `phi=0.8` is their phase difference in radians. To identify
the global ground state, compare the relevant even and odd sectors.

The bath and many-body settings are:

- `pairs=1`: one mirrored pair of normal-state bath energies, at $+\xi$ and
  $-\xi$, per lead. Increasing `pairs` refines the bath discretization; it does
  not count leads or spin directions.
- `compress=False`: retain all bath modes. In this symmetric model, some linear
  combinations of degenerate bath QP modes have tunneling amplitudes to the
  impurity that cancel exactly. These "dark" combinations decouple from the
  impurity, although occupying them still costs excitation energy.
  `compress=True` retains the coupled ("bright") combinations and freezes the
  discarded dark modes in their vacuum, omitting states with dark-mode
  excitations. Use `False` when the full finite-bath spectrum is needed.
- `parity=1`: select odd total fermion parity; `0` selects even parity. This
  includes the impurity and both reservoirs, not just the impurity.
- `twice_sz=1`: select total spin projection $S_z=+1/2$, in units where
  $\hbar=1$. This also includes both impurity and reservoir spins.
- `cutoff=2`: allow at most two bath QPs in total across both leads and spins.
  All impurity configurations remain available, subject to the selected total
  sector. `cutoff=None` removes this bath-QP truncation.

The printed output is approximately:

```text
[-0.2950596] [[0.89270561 0.10104019 0.00625420]]
```

`qp.energies` contains the computed many-body energies in the selected sector,
not excitation gaps. By default, `solve` returns only the lowest state.
`qp.qp_weights[i, n]` is the probability of exactly `n` bath QPs in eigenstate
`i`, obtained by summing squared amplitudes of the corresponding configurations;
impurity occupation is not counted. Here the single row gives probabilities
of approximately **89.27%, 10.10%, and 0.63%** for zero, one, and two bath QPs.
Each row sums to one. Odd total parity does not require odd bath-QP number,
because the impurity also contributes to parity.

The assertions check the reference energy, the eigenvalue-equation residual
`qp.residuals`, and normalization of the QP weights. The residual measures
eigensolver accuracy within the truncated space; neither a tiny residual nor
normalized weights establishes convergence with QP cutoff or bath size.

To continue with DMRG, the `[dmrg]` extra must be installed. If you chose the base
installation, first run `python -m pip install '.[dmrg]'` from the main source
directory, then continue in the same Python session:

```python
import numpy as np

mps = solve(h, sector=sector, backend="dmrg",
            options={"eigenpairs": 2, "chi_max": 32, "require_convergence": True})
exact = solve(h, cutoff=None, sector=sector, options={"eigenpairs": 2})
np.testing.assert_allclose(mps.energies, exact.energies, rtol=0, atol=1e-9)
print(mps.energies, mps.residuals)
print(mps.matrix_elements(h.observables["impurity_spin_z"]))
```

Both methods use the same finite Hamiltonian and energy zero. The first
QP calculation, `qp`, uses cutoff two; `exact` includes every configuration
allowed by the selected symmetry sector. DMRG includes the same modes and
approximates the state through its MPS bond dimension. Eigensolver accuracy,
QP/MPS convergence, and bath
convergence must be assessed separately.

To save a result and its wavefunctions, then load them for further analysis:

```python
from qdjj_solver import save_result, load_result

save_result(qp, "results/state.json", save_states=True)
restored = load_result("results/state.json")
```

Python `solve` includes all configurations of the specified finite bath by
default; the `qdjj-solver` terminal command defaults to QP cutoff two.
Set `compress=False` explicitly when constructing a full reference
model in Python. See the
[defaults table](https://github.com/rokzitko/qdjj_solver/blob/main/docs/interfaces.md#defaults)
and [convergence guidance](https://github.com/rokzitko/qdjj_solver/blob/main/docs/numerics.md)
before increasing problem sizes.

The [`first_calculation/`](https://github.com/rokzitko/qdjj_solver/blob/main/first_calculation/README.md)
directory provides JSON inputs and reference outputs for these calculations.
Run `python -B first_calculation/run.py` for QP with the base installation.
With the `[dmrg]` extra installed, add `--all` for all three
calculations (QP-expansion, DMRG, exact diagonalisation).

## Command line

To specify a calculation in an input file, save the following text as
`model.json` in your working directory. This example requires only the base
installation:

```json
{
  "format_version": 1,
  "backend": "qp",
  "model": {"kind": "reference", "u": 1.6, "gamma": 0.3, "phi": 0.8, "compress": false},
  "bath": {"kind": "cosh-grid", "pairs": 1, "bandwidth": 5.0},
  "cutoff": 2,
  "sector": {"parity": 1, "twice_sz": 1},
  "solver": {"eigenpairs": 2}
}
```

```sh
qdjj-solver solve model.json --output results/qp.json --save-states
```

The command creates the output directory and prints a lowest energy near
`-0.2950595967`. It saves the energies and calculation settings in JSON, and the
wavefunction coefficients in a companion NumPy file (NPZ). Load the result with
`load_result("results/qp.json")`.

### Repository examples and DMRG continuation

The following commands use input files supplied with the source code; run them
from the main source directory. The DMRG solve and `--resume` commands require
the `[dmrg]` extra: install it first with `python -m pip install '.[dmrg]'`.

```sh
qdjj-solver solve examples/qp.json --output results/qp.json --save-states
qdjj-solver solve examples/dmrg.json --output results/dmrg.json --save-states
qdjj-solver solve examples/dmrg_refined.json --resume results/dmrg.json --output results/refined.json --save-states
```

`dmrg_refined.json` raises `chi_max` from 32 to 64, permits more sweeps, and tightens
energy and residual tolerances. `--resume` starts a new optimization from the
saved, completed MPS states using the new settings. This tiny example normally
already converges at bond dimension 32; the same procedure is useful for larger models.

`qdjj-solver` saves both methods' results in a common JSON format. Wavefunctions
are saved separately as NPZ files for QP calculations and TeNPy HDF5 files for
DMRG. Keep the JSON and its named wavefunction file together when copying a
calculation.

## Documentation

### Getting Started

- [First calculation: inputs and reference outputs](https://github.com/rokzitko/qdjj_solver/blob/main/first_calculation/README.md)
- [Runnable examples](https://github.com/rokzitko/qdjj_solver/blob/main/examples/README.md)
- [Worked parity, excitation, and current calculation](https://github.com/rokzitko/qdjj_solver/blob/main/docs/numerics.md#worked-parity-and-current-example)

### Reference

- [Models, units, sectors, and observables](https://github.com/rokzitko/qdjj_solver/blob/main/docs/models.md)
- [Python functions, input files, saved results, and troubleshooting](https://github.com/rokzitko/qdjj_solver/blob/main/docs/interfaces.md)
- [Validation and development](https://github.com/rokzitko/qdjj_solver/blob/main/docs/validation.md)
- [Release guide for maintainers and developers](https://github.com/rokzitko/qdjj_solver/blob/main/docs/releases.md)
- [Citation and method references](https://github.com/rokzitko/qdjj_solver/blob/main/docs/references.md)

### Convergence

- [Bath representations: parameters, convergence, and precision/time benchmarks](https://github.com/rokzitko/qdjj_solver/blob/main/docs/bath_representations.md)
- [Numerical methods and convergence](https://github.com/rokzitko/qdjj_solver/blob/main/docs/numerics.md)
- [QP-cutoff convergence: energy, local moment, charge probabilities, and dot–lead spin correlation](https://github.com/rokzitko/qdjj_solver/blob/main/docs/qp_convergence.md)

### Research Results

- [Large-hybridization study: 2QP, 4QP, 6QP and DMRG through Gamma/Delta = 10](https://github.com/rokzitko/qdjj_solver/blob/main/large_Gamma/README.md)
- [Research applications and literature reproductions](https://github.com/rokzitko/qdjj_solver/blob/main/applications/README.md)
- [NRG comparison results](https://github.com/rokzitko/qdjj_solver/blob/main/NRG_comparisons/README.md) and [test1 results](https://github.com/rokzitko/qdjj_solver/blob/main/NRG_comparisons/test1/README.md)

## Support and contributing

Before opening a [GitHub issue](https://github.com/rokzitko/qdjj_solver/issues),
consult the [troubleshooting and issue-report checklist](https://github.com/rokzitko/qdjj_solver/blob/main/docs/interfaces.md#troubleshooting-and-exit-status).
See [CONTRIBUTING.md](https://github.com/rokzitko/qdjj_solver/blob/main/CONTRIBUTING.md)
for development setup and the test, lint, and packaging commands used by CI.
Reinstall after changing the native C++ implementation or build configuration.
