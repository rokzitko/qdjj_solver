# Small example calculations

Install from source following the root README's
[installation prerequisites](https://github.com/rokzitko/qdjj_solver/blob/main/README.md#install),
including the native C++ toolchain and Python development headers. Run all
commands below from the main source directory.

## Base-only examples

These commands need only the base installation (`python -m pip install .`);
TeNPy and HDF5 are not required:

```sh
python -B examples/quickstart.py
python -B examples/sector_convergence.py
python -B examples/qp_convergence.py
qdjj-solver solve examples/qp.json --output results/qp.json --save-states
```

`quickstart.py` needs only the base installation and checks the README's
finite-model energy, residual, and normalization. `sector_convergence.py` is also
QP-only; it compares both parities, refines the QP cutoff, and computes an
excitation energy and current. See its [worked interpretation](https://github.com/rokzitko/qdjj_solver/blob/main/docs/numerics.md#worked-parity-and-current-example).
`qp_convergence.py` follows a spin-up doublet through a QP-cutoff ladder to full
finite-bath ED, reporting energy, local spin, dot charge probabilities, and a
direct dot–bath spin correlation. It needs only the base installation and checks
the SU(2) doublet identity at every cutoff. The
[cutoff convergence guide](https://github.com/rokzitko/qdjj_solver/blob/main/docs/qp_convergence.md) adds a larger reference study.

## DMRG-required examples

Before running these commands, install the `[dmrg]` extra with
`python -m pip install '.[dmrg]'`. It supplies TeNPy and HDF5 support and is
required for both the DMRG commands and `multi_orbital.py`:

```sh
qdjj-solver solve examples/dmrg.json --output results/dmrg.json --save-states
qdjj-solver solve examples/dmrg_refined.json --resume results/dmrg.json --output results/refined.json --save-states
python -B examples/multi_orbital.py
```

`qp.json` and `dmrg.json` use identical physical coordinates. The QP
example imposes a finite cutoff; change it to `null` for unrestricted finite
space. The multiorbital example compares unrestricted ED and DMRG for a small
spin-mixing model. These examples show how to set up and check a calculation;
their small baths serve as finite-model illustrations.

`dmrg_refined.json` uses the same physical model and MPS site ordering as
`dmrg.json`, raises the maximum bond dimension to 64, and tightens the numerical
tolerances. The `--resume` command
initializes a new optimization from the first run's saved MPS states. Both
commands save wavefunctions so either completed result can be loaded again.

These files are supplied with the source code. To run outside the source tree,
the [main README](https://github.com/rokzitko/qdjj_solver/blob/main/README.md#command-line)
provides a complete JSON input that can be saved in any working directory.
