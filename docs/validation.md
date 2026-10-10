# Validation and development

## Local checks

```sh
python -B tools/check_publication.py
python -m pip install -e '.[dmrg,dev]'
python -B -m pytest --cov --cov-report=term-missing --cov-report=json
python tools/check_coverage.py reports/python-coverage.json
python -m ruff check src tests tools examples applications first_calculation large_Gamma NRG_comparisons
python -m build
python tools/check_distribution.py dist/*.tar.gz dist/*.whl
python -B tools/check_publication.py dist/*.tar.gz dist/*.whl
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

## Publication boundary

`pyproject.toml` and `tools/check_distribution.py` contain matching explicit NRG
asset allowlists, checked by `tests/test_ci_tools.py`. Test2's asset list is
unchanged. Test3 includes its standard-library generator and
`reference/provenance.json`. Test4 ships exactly 22 reviewed files: model sources,
generic numerical/reproduction helpers, two physical/bath inputs, documentation,
and the five top-level scientific outputs. `FINAL_REPORT.md` is only a pointer
to the current report. Deleted operational controllers, old run profiles and
intermediate output subdirectories remain negative packaging fixtures even
though they no longer exist in the public tree.

Public Test2/3/4 JSON uses an explicit scientific-field whitelist, not a recursive
dump of execution records. Retain model and bath definitions, numerical settings,
state observables, convergence qualifications, failed or unavailable comparisons,
and appropriately scoped checksums. Test4's normalized snapshot retains 157
records in approximately 4.65 MiB, including the accumulated negative scientific
findings. Repeated scientific definitions are referenced rather than duplicated.
Private operational archives, controllers and execution identities are not part
of this schema. Original-artifact hashes identify private originals; they do not
claim that those originals were revalidated from the public snapshot. Removing
historical public copies does not remove or replace privately retained raw runs.

Regenerate the public tables and reports into fresh directories using only the
standard library, without solver imports, raw runs, network access or private
archives:

```sh
python -B -S NRG_comparisons/test2/generate.py --tables-only --output-dir results/test2-presentation
python -B -S NRG_comparisons/test3/generate.py --output-dir results/test3-presentation
python -B -S NRG_comparisons/test4/generate.py --tables-only --output-dir results/test4-presentation
```

These commands replay published data, not expensive numerical calculations or a
new continuum-accuracy certification. Test2's optional PNG output still requires
the plotting extra. Test1's assets and existing dataset policy are unchanged;
Test1 is exempt from the new semantic-metadata and dataset-size policies.

The read-only, offline `tools/check_publication.py` complements the filename
boundary with semantic privacy checks. With no arguments it checks tracked and
unignored untracked files under Test2/3/4, including their reviewed manual code;
it does not traverse ignored raw runs. Explicit source paths work without Git.
The gate rejects private paths, site/account literals, operational JSON metadata,
credential patterns and oversized public data. Default limits are **8 MiB per
JSON file** and **16 MiB per case**. Schema validation remains the generators' and
case tests' responsibility; a successful privacy scan alone is not scientific
validation.

Archive targets are inspected without extraction or ignore exclusions. Private
members are rejected, and all archive text is checked for credential patterns,
including files outside Test2/3/4. Semantic metadata and dataset-size checks
remain scoped to Test2/3/4. Binary files are counted as skipped, not certified
secret-free. Findings produce a failing exit status with redacted diagnostics;
no external SaaS secret scanner or upload is involved.

The source-root allowlist is deliberately broader outside `NRG_comparisons`, so
filename acceptance alone is not a privacy guarantee. Keep private archives,
worktrees, Git bundles, virtual environments, credentials and raw `runs/` out of
public assets. Both the distribution boundary and publication gate must pass
on the actual archives before upload, not just on the checkout.

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

Every test job runs the default source privacy gate before installation, lint,
examples or tests can print dataset contents. Diagnostic report uploads are
disabled if that gate fails. Built distributions pass archive privacy checks
after the distribution-boundary check and before artifact upload. Release jobs
also check their source checkout, recheck the downloaded sdist before unpacking,
and check repaired wheels before upload. These are mandatory local checks in CI,
not external scanning services.

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
