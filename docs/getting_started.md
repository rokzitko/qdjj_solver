# Installation and a first checked calculation

This chapter explains how to install the package from source and calculate a quantum-impurity
state. It assumes basic Python, including imports, function calls and dictionaries.
The numerical methods are introduced as they are used. The physics
manual explains the Hamiltonian and subgap states; this handbook explains how to
set up and reproduce the corresponding numerical calculation.

## Install the package

Use Python 3.11 or later. Building from source needs a C++17 compiler and the
development headers for that Python interpreter. On Debian/Ubuntu these normally
come from `build-essential` and the matching `python3-dev`; macOS uses the Xcode
Command Line Tools; Windows uses Visual Studio C++ Build Tools and the Windows
SDK. Install the compiler and headers first; pip supplies CMake and pybind11.

```sh
git clone https://github.com/rokzitko/qdjj_solver.git
cd qdjj_solver
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or
`.venv\Scripts\Activate.ps1` in Windows PowerShell. Choose one installation:

```sh
python -m pip install .
# Or include the optional DMRG backend and plotting:
python -m pip install '.[dmrg,plots]'
```

The base installation supports the model-building API and QP solver. The `dmrg`
extra adds TeNPy and support for saving MPS states in HDF5 files; `plots` adds Matplotlib. Development and
documentation dependencies are only needed for their corresponding tasks.
Run the following command in that same environment to check which package is used:

```sh
python -c "import qdjj_solver; print(qdjj_solver.__version__, qdjj_solver.__file__)"
```

Commands below that begin with `examples/` are run from the source root. The
example scripts resolve their supplied inputs relative to their own locations,
so an absolute script path also works from a separate working directory.

## Calculate one state

The following code is generated from the complete, independently runnable
[`examples/quickstart.py`](../examples/quickstart.py). Copy it into a Python
session, or run `python -B examples/quickstart.py`.

<!-- BEGIN generated: quickstart (tools/gen_docs.py) -->
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
<!-- END generated: quickstart -->

The output contains one energy followed by the probabilities of zero, one and
two bath quasiparticles. They are approximately 89.27%, 10.10% and 0.63% here.
The script checks the energy against its finite-model reference, the residual
of the eigenvalue equation, and normalization of those probabilities.

An **object** stores related information: `bath` holds a discretized reservoir,
`h` holds the finite Hamiltonian, `sector` identifies allowed quantum numbers,
and `qp` holds the numerical result. A **backend** is the numerical method used
to solve that same Hamiltonian. The top-level `solve` selects one with `backend`.

| Input | Meaning in this calculation |
|---|---|
| `pairs=1` | Two signed normal levels per lead, each with two spin modes |
| `delta=1` | Energy unit: the superconducting gap |
| `bandwidth=5` | Normal-state half-bandwidth |
| `gamma=0.3` | Total hybridization, shared equally by the two leads |
| `compress=False` | Retain the full declared finite bath, including dark modes |
| `Sector(1, 1)` | Odd total parity and total spin projection +1/2 |
| `cutoff=2` | At most two bath QPs; all impurity configurations are retained |

The reported energy is the many-body energy of a state in the selected sector.
Excitation energies are differences between state energies, and identifying the global
ground state requires comparison with the other relevant sectors. The projected residual
measures eigensolver accuracy within the chosen space; cutoff and bath refinement assess
the accuracy of the approximations defining that space.

## Run the same case through an input file

Tutorial 00 supplies complete inputs for QP, unrestricted same-bath ED and DMRG:

```sh
qdjj-solver solve examples/00_first_calculation/input/qp.json --output results/first.json --save-states
```

This creates a JSON record and a companion NPZ file containing the wavefunctions.
Inspect the JSON filename printed by the command and the saved diagnostics.
With DMRG, the convergence flags identify whether the returned states satisfy the
finite-problem checks; its default also permits diagnostic output from unconverged runs.

Python and CLI defaults differ. In particular, Python `solve` defaults to
`cutoff=None`, while the QP CLI defaults to two. Always set consequential numerical
and representation choices explicitly in a scientific calculation. The full
[defaults table](interfaces.md#defaults) lists these differences.

## Compare methods and keep the result

With `.[dmrg]` installed, run the full first-calculation tutorial:

```sh
python -B examples/00_first_calculation/run.py --profile full --output runs/first-comparison
```

This compares the same finite Hamiltonian using QP truncation, DMRG and full
finite-space ED. Agreement with ED validates the solution of that finite problem.
Assess its accuracy as a continuum approximation through separate bath refinement.

Use the [saved-result workflow](workflows.md) to distinguish scalar-only analysis
from reloading wavefunctions. The latter requires keeping the JSON and its named
companion file together. A **checkpoint** here is a completed saved MPS state
used as the initial wavefunction for a new optimization.

## Choose your next route

- To calculate subgap states and current, read [model construction](models.md),
  [running the solver](interfaces.md), then [static calculations](static_calculations.md).
- To establish accuracy, read [reservoir refinement](reservoir_workflow.md),
  [QP convergence](qp_workflow.md) or [DMRG refinement](dmrg.md).
- To calculate driven spectra, first obtain and interpret the static states,
  then follow [Floquet calculations](floquet.md) and [spectroscopy](spectroscopy.md).
- To look up a parameter, use the [constructor reference](reference/models.md),
  [solver options](reference/options.md) and [input-format reference](reference/configuration.md).

The small installation test above and the reference junction used in the next
workflow are deliberately different cases. Compare expected energies using
matched Hamiltonians and numerical settings.
