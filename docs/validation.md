# Validation and development

## Local checks

```sh
python -m pip install -e '.[dmrg,dev]'
python -B -m pytest --cov --cov-report=term-missing --cov-report=json
python tools/check_coverage.py reports/python-coverage.json
python -m ruff check src tests tools examples applications first_calculation large_Gamma NRG_comparisons
python -m build
python tools/check_distribution.py dist/*.tar.gz dist/*.whl
python tools/check_wheel.py dist/*.whl --dmrg
```

The regular suite excludes `slow` tests. Run the bounded extended checks with
`python -B -m pytest -m slow`. Numerical tests use fixed seeds and explicit
tolerances; small dense electron-space matrices and analytic spectra provide
independent references. Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` to
reproduce CI's thread limits.

The wheel check installs the exact artifact and resolves its dependencies in a
fresh virtual environment outside the checkout. It runs `pip check`, verifies
package/extension origins, exercises the `qdjj-solver` command, and checks
numerical results and serialization. It also blocks TeNPy imports during its
QP portion.
The archive review checks the software-only source boundary and required members.
For the [NRG comparison](../NRG_comparisons/README.md), the source archive includes
the reviewed results and an explicit allowlist of reproduction scripts and
parameter inputs, including the seven final finite baths. Raw NRG outputs,
checkpoints and historical run directories remain excluded. These study assets
are not included in wheels.

The NRG adapter/control and reproduction tests are bounded and offline: literal
input decks, synthetic HDF5, small electron-ED/QP problems, task/input consistency,
and fresh-directory preparation. Optional plotting tests regenerate presentation
files in temporary directories without changing the approved output. Follow the
[case reproduction instructions](../NRG_comparisons/test1/README.md#reproduction)
for expensive QP/DMRG and external NRG runs. Those require Linux and, for NRG,
NRG Ljubljana, Mathematica and `h5py`; routine tests do not launch them or certify
the published precision on a new installation.

## CI and coverage policy

`.github/workflows/tests.yml` has five routine jobs:

| Job | Checks |
|---|---|
| Linux, Python 3.13 | Both backends, lint, Python coverage, sdist-to-wheel build, fresh-wheel smoke test |
| macOS, Python 3.14 | Native/QP, persistence, command-line checks, platform wheel |
| Windows, Python 3.11 | Native/QP, persistence, command-line checks, platform wheel |
| Linux, Python 3.13 | Focused native ASan/UBSan regressions |
| Linux, Python 3.11 | Both backends with the runtime lower bounds in `tools/minimum-runtime.txt` |

The Monday schedule and manual `workflow_dispatch` also run the stress checks
and a native-coverage diagnostic job. Jobs have timeouts, single-thread numerical
defaults, and 14-day diagnostic artifact retention. Superseded PR runs are cancelled.
Routine wheel/sdist artifacts are retained for 14 days as well. The release
workflow calls the same checks with stress and native coverage enabled; see
[the release procedure](releases.md).

The lower-bound job constrains all five numerical/runtime dependencies together,
runs `pip check`, requires the installed-package smoke check (including DMRG and
HDF5 persistence), and executes the regular suite. Missing optional backends fail
the smoke check rather than silently skipping their validation.
To reproduce it in a fresh Python 3.11 environment:

```sh
python -m pip install -c tools/minimum-runtime.txt '.[dmrg,test]'
python -m pip check
python -B tools/check_wheel.py --installed --dmrg
python -B -m pytest
```

The required coverage floor is **95% Python executable lines**, including the
canonical and compatibility packages. `coverage.py` measures subprocesses and
unimported modules, with installed/source paths mapped for reporting. The small
`tools/check_coverage.py` gate uses line counts rather than the combined
line/branch percentage shown by `coverage report`. Python branch coverage and
C++ coverage are reported for review. Test code and third-party implementations
are outside the product coverage denominator.

Add `--cov-report=html --cov-report=xml` for browsable and machine-readable
reports under `reports/`. Coverage is evidence about executed code; the numerical
assertions and independent oracles determine what was actually verified.

## Native instrumentation

Use a fresh build directory for each instrumented run. On Linux with GCC:

```sh
python -m pip install 'scikit-build-core>=0.11' 'pybind11>=3.0' cmake ninja 'gcovr>=8.4'
python -m pip install --no-build-isolation -e '.[test]' \
  -Cbuild-dir=build/native-coverage -Ccmake.build-type=Debug \
  -Ccmake.define.QDJJ_COVERAGE=ON
python -B -m pytest tests/test_native_core.py tests/test_qp.py --junitxml=reports/native-tests.xml
python -m gcovr --root . --filter 'cpp/core.cpp' \
  --gcov-object-directory build/native-coverage \
  --exclude-throw-branches --exclude-unreachable-branches \
  --html-details reports/native.html --print-summary
```

This build keeps the pybind11 headers available for gcov's source lookup.
On macOS with Apple Clang, add `--gcov-executable 'xcrun llvm-cov gcov'` to
the reporting command. Compiler-generated exception/unreachable arcs are
filtered; native branch percentages remain compiler-dependent diagnostics.

`QDJJ_SANITIZERS=ON` enables address and undefined-behavior instrumentation in
a separate build. The workflow shows the matching GCC runtime preload needed
when importing an instrumented extension into ordinary Python. Address and UB
failures are fatal. Process-wide leak detection is disabled because CPython and
third-party libraries are outside the extension's allocation-lifetime contract.
Reinstall normally after local instrumentation:

```sh
python -m pip install -e '.[dmrg,dev]'
```

## Independent physics checks

- Independent tensor-product physical-electron ED versus QP and DMRG, including
  complex reservoir contacts, spin mixing, and energy-zero conventions.
- Multiorbital interactions, exchange, pair hopping, and orbital mixing.
- Random fermion strings versus tensor products; projection after multiplication.
- Fermionic exchange signs for large mode counts and exact sector dimensions.
- Orthogonal excited states with positive energies, additional conserved subspaces, arbitrary
  mode permutations, grouping choices, eta and particle-number sectors.
- Matrix elements between parity/spin sectors, overlap matrices, full residuals,
  saved MPS wavefunctions, and continued optimization.
- Long finite MPS representations, distinguishing mode ordering along the MPS
  from physical geometry.
- Quadrature normalization, hybridization functions, published surrogate residues,
  finite-difference physical derivatives, QP cutoff monotonicity, and eigensolver agreement.
- Consistency of older and current Python functions and saved results, plus
  `qdjj-solver` calculations run outside the source directory.
- Literature applications: archived finite-bath baselines, perturbative Knight
  shift, Zeeman/spin consistency, double-dot current derivatives and GAL phase
  boundaries, and direct three-lead versus effective two-lead energies/currents.

These checks establish numerical accuracy for the tested finite Hamiltonians.
Continuum accuracy of a new
physical calculation requires its own bath and solver convergence study as
described in [numerics](numerics.md). DMRG checks are skipped when TeNPy is absent;
the model and QP checks remain usable.

Additional tests cover invalid inputs to the C++ routines, storage of basis
states, Schur solutions near eliminated-block poles, complex transition matrix
elements, and interrupted saves of wavefunctions. Earlier calculation and
file-reading routines receive representative physical and usage checks.

Complex degenerate and near-degenerate QP multiplets are checked against dense
one-particle spectra, including calculations that automatically select sparse
diagonalization. Input-file tests also check invalid JSON, misspelled settings,
and mismatched result/wavefunction files.

## Checking the research applications

The [application guide](../applications/README.md#numerical-checks-and-further-studies)
lists commands for small checks of all ten studies. Reference-value checks reuse
the stored bath nodes, so changes in the least-squares optimizer cannot change
the finite Hamiltonian being compared. Other checks fit new baths, run the
calculation scripts, and save JSON/CSV results in temporary directories.
The extended checks selected with `-m slow` include small QP/DMRG comparisons.
These tests need neither plotting nor network access.

The `paper` and `convergence` profiles produce scientific reports and can take
minutes or hours. Run them explicitly for a literature reproduction or a
convergence study. The archived outputs document the achieved literature
agreement and convergence; the small finite-bath checks alone do not establish
continuum accuracy. The source archive includes the application scripts, inputs,
reference data, and compact results. The precompiled package provides the solver
used to run these calculations.

The [large-hybridization study](../large_Gamma/README.md) has bounded independent
quadratic/ED and optional DMRG checks in `tests/test_large_gamma.py`. Its full
multi-day scan is an explicit command, with compact reviewed outputs and
working records under ignored `runs/large_Gamma/`.
