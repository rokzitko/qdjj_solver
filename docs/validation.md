# Validation and development

## Local checks

```sh
python -B tools/check_publication.py
python -m pip install -e '.[dmrg,dev]'
python -B -m pytest --cov --cov-report=term-missing --cov-report=json
python tools/check_coverage.py reports/python-coverage.json
python -m ruff check src tests tools examples applications large_Gamma NRG_comparisons docs/conf.py docs/_ext docs/manual/scripts
python -m build
python tools/check_distribution.py dist/*.tar.gz dist/*.whl
python -B tools/check_publication.py dist/*.tar.gz dist/*.whl
python tools/check_wheel.py dist/*.whl --dmrg
```

The wheel checker accepts **one wheel**. If `dist/` contains several, replace `dist/*.whl`
in the last command with one wheel compatible with your interpreter and platform; keep
`--dmrg` to check both backends and HDF5 persistence.

Documentation changes, and code changes that affect documented behaviour, also need the
documentation checks (install `docs`, `docs-pdf` and `plots` for the PDF checks);
see [documentation](documentation.md#checks) and [PDF build requirements](user_manual.md#build-locally):

```sh
python tools/gen_docs.py --check
python tools/check_docs.py
python -B -m pytest -m "docs or tutorial"
python -m sphinx -W --keep-going -c docs -b html . build/docs
python -B tools/build_manuals.py --manual all
```

The regular suite excludes the `slow`, `docs` and `tutorial` tests. Run the
extended checks with `python -B -m pytest -m slow`. Numerical tests use fixed seeds and
explicit tolerances; small dense electron-space matrices and analytic spectra provide
independent references. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`,
`MKL_NUM_THREADS=1` and, on macOS, `VECLIB_MAXIMUM_THREADS=1` to reproduce the single
numerical thread of the CI jobs (continuous integration: the checks that GitHub runs for
every pull request).

The wheel check installs the built wheel and its dependencies in a new virtual
environment outside the repository. It runs `pip check`, verifies that the package and
its compiled extension are imported from that installation, runs the `qdjj-solver`
command, and checks numerical results and saved files. TeNPy imports are blocked during
its QP part, which must work without TeNPy.
`tools/check_distribution.py` checks that the source archive contains the required files
and only files from its explicit lists.
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

The tests of NRG input preparation and reproduction use small, offline cases: literal
input decks, synthetic HDF5, small electron-ED/QP problems, task/input consistency,
and fresh-directory preparation. Optional plotting tests regenerate presentation
files in temporary directories without changing the approved output. Follow the
[case reproduction instructions](../NRG_comparisons/test1/README.md#reproduction)
for expensive QP/DMRG and external NRG runs. Those require Linux and, for NRG,
NRG Ljubljana, Mathematica and `h5py`; routine tests do not launch them or establish
the published precision on a new installation.

## CI and coverage policy

`.github/workflows/tests.yml` has six routine jobs:

| Job | Checks |
|---|---|
| Linux, Python 3.13 | Both backends, `ruff`, Python coverage, source distribution and wheel build, test of the freshly installed wheel |
| Linux, Python 3.13, documentation | `tools/check_docs.py`, the `docs` and `tutorial` tests, strict HTML and user/reference PDF builds |
| macOS, Python 3.14 | Native/QP, persistence, command-line checks, platform wheel |
| Windows, Python 3.11 | Native/QP, persistence, command-line checks, platform wheel |
| Linux, Python 3.13 | Focused native ASan/UBSan regressions |
| Linux, Python 3.11 | Both backends with the runtime lower bounds in `tools/minimum-runtime.txt` |

The weekly run (Mondays) and manual runs (`workflow_dispatch`) add the stress checks
(`-m slow`), a native-coverage diagnostic job and the extended documentation job: the
two PDF manuals, physics listings, the full tutorial profiles and the recomputed manual
data. Jobs have time limits and use one numerical thread. Their reports and wheels are kept
for 14 days, the PDF artifacts for 30 days. A new push to a pull request cancels its
running checks. The release workflow runs the same checks, including the weekly ones; see
[the release procedure](releases.md).

Every test job runs the default source privacy gate before installation, lint,
examples or tests can print dataset contents. Diagnostic report uploads are
disabled if that gate fails. Built distributions pass archive privacy checks
after the distribution-boundary check and before artifact upload. Release jobs
also check their source checkout, recheck the downloaded sdist before unpacking,
and check repaired wheels before upload. These are mandatory local checks in CI,
not external scanning services.

The lower-bound job constrains all five numerical/runtime dependencies together,
runs `pip check`, requires the installed-package check (including DMRG and
HDF5 persistence), and executes the regular suite. Missing optional backends fail
the installed-package check rather than silently skipping their validation.
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
are excluded from the package's coverage totals.

Add `--cov-report=html --cov-report=xml` for browsable and machine-readable
reports under `reports/`. Coverage is evidence about executed code; the numerical
assertions and independent reference calculations determine what was actually verified.

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

`QDJJ_SANITIZERS=ON` enables address and undefined-behaviour instrumentation in
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
- Floquet quasienergies (`tests/test_floquet.py`): Jacobi–Anger harmonics of the phase
  drive; band-limited linear drives; the coordinate and sector checks; the undriven limit;
  Bessel weights of a decoupled dot with both scalar gauges; the exactly solvable circularly
  driven two-level system (Sambe, propagator, RWA and the analytic van Vleck operators);
  van Vleck error scaling with $\Omega$; Sambe versus propagator in one basis; static-basis
  convergence to the Fock result; the rigorous truncation bound in $M$; the phase–voltage
  gauge identity; the RWA splitting at the resonance between $S$ and $S^*$; the eta-sector
  decomposition; and
  propagation of the exact $H(t)$ of a phase- and detuning-driven junction in the physical
  electron space, with a different Magnus integrator.
- Probe matrix elements and spectra (`tests/test_floquet_transitions.py`): static elements
  against an independent electron-space calculation; Floquet elements against independent
  time propagation (with a wrong-sign control); the sum rule, Hermiticity and replica invariance; Tien–Gordon
  sidebands across parity sectors in both scalar gauges; analytic lines of the circularly
  driven two-level system; the undriven limit with basis-invariant merging of degenerate
  states; agreement of the Sambe, propagator and static-basis methods; and the adiabatic
  generalized Tien–Gordon limit of a gate-driven junction.
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
reference data, and compact results; the wheels contain only the solver packages
used to run these calculations.

The [large-hybridization study](../large_Gamma/README.md) has small independent
quadratic/ED and optional DMRG checks in `tests/test_large_gamma.py`. Its full
multi-day scan is an explicit command, with compact reviewed outputs and
working records under ignored `runs/large_Gamma/`.
