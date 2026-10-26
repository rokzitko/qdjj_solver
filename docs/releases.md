# Release artifacts and procedure

This guide describes release preparation for **maintainers and developers**.
For installation from source, follow the
[root README](../README.md#install).
The PyPI publication steps below are for a future release; they do not imply
current PyPI availability.

## Versioned artifacts

`.github/workflows/release.yml` runs on `v*` tags and manual dispatch. A manual
run rehearses the full process and retains artifacts without publishing a
GitHub release, including when dispatched against a tag. Only a tag **push**
attaches the tested sdist, wheels and both PDF manuals to a GitHub release after
all checks succeed.

The workflow:

1. Runs the shared CI suite, including minimum runtime dependencies, sanitizers,
   stress tests, native coverage and the extended documentation checks.
2. Checks agreement between `pyproject.toml`, `__version__`, and `CITATION.cff`;
   a release tag must be exactly `v<version>`.
3. Builds and inspects one source distribution.
4. Builds wheels from that sdist with cibuildwheel for CPython 3.11–3.14:
   Linux x86-64/manylinux 2.28, macOS x86-64 and arm64, and Windows AMD64.
   macOS extension builds target deployment version 11.0; runtime dependency
   wheels can have their own OS requirements.
5. Exercises each installed wheel with both numerical backends, state
   serialization, optional-import isolation, and console commands. Cibuildwheel
   supplies the isolated test environment; `tools/check_wheel.py --installed
   --dmrg` checks the installed artifact there.
6. Runs archive-boundary and strict package-metadata checks and retains release
   artifacts for 30 days. Tagged releases attach those same files permanently.
7. Installs the exact source archive in a documentation job, runs the small
   handbook examples, and builds both manuals with the reviewed documentation
   toolchain. The PDFs and `qdjj-solver-<version>-manuals.json` record their
   source commit/archive, software identity, build tools, page counts and checksums.

Manual filenames are `qdjj-solver-<version>-physics-manual.pdf` and
`qdjj-solver-<version>-users-guide-reference.pdf`. Publication requires both
validated artifacts and their manifest. Build logs are retained separately on
failure. The [PDF guide](user_manual.md) documents local builds and system tools.

## Preparing a release

Before making the repository public, perform a **clean-history publication step**
from the reviewed source tree. The private development history includes local
machine paths in earlier research archives; deleting those files does not remove
them from older commits. Keep the development history privately, and validate
the clean public candidate and its freshly built distributions before publishing.
This is a pending publication step, not an instruction to rewrite history during
routine development.

Update the three version declarations together, refresh the release date and
the citations in [references](references.md#software-citation) and
[`qdjj_solver.bib`](qdjj_solver.bib), then check the declarations locally:

```sh
python -B tools/check_release.py --tag v0.3.0
```

Use the all-versions Zenodo DOI until the new version-specific archive is
available, then update the DOI in `CITATION.cff` and both citation documents.

Run a manual release-artifact workflow on the candidate commit and inspect all
job results. Confirm the public installation instructions and runnable examples
match the candidate. Create the matching tag on the verified commit when ready.
Verify that the documentation and artifact links of the release work without
authentication.

For example, to rehearse the current `main` branch:

```sh
gh workflow run release.yml --ref main
```

## Future PyPI publication

For a future PyPI release, use the GitHub release's already tested sdist and
wheels as the final distributable files. Install Twine first, then download
only that release's `.tar.gz` and `.whl` assets into a new, otherwise empty
`dist/pypi-upload/` directory. Run these commands from the repository root;
publication requires authorized PyPI credentials:

```sh
python -m pip install twine
# Download the tested GitHub release assets into dist/pypi-upload/ before continuing.
python -m twine check --strict dist/pypi-upload/*.tar.gz dist/pypi-upload/*.whl
python -m twine upload dist/pypi-upload/*.tar.gz dist/pypi-upload/*.whl
```

Upload exactly the downloaded, checked files; do not rebuild them for PyPI.
The release workflow itself does not publish to PyPI or require PyPI credentials.
Avoid rebuilding different files under an already published version.

Before announcing a future PyPI release, also verify that the advertised PyPI
installation resolves the published version. This additional check does not
apply to the current source-build installation instructions.

## Release notes

### Unreleased

- Periodically driven junctions: Floquet states from the Sambe matrix (dense or
  shift-invert), the one-period propagator, the van Vleck expansion and the rotating-wave
  approximation, with residual bounds and convergence diagnostics.
- Matrix elements of probes between Floquet states, and between static and Floquet states;
  spectral lines with golden-rule strengths, sum rules and broadened spectra.
- The physics manual (LaTeX), this user and developer manual (Sphinx), and
  runnable tutorials with reference results, kept consistent by generated reference pages
  and documentation checks in CI.
- A separate User's Guide and Reference Manual PDF with task-oriented workflows,
  generated constructor/configuration/record references, and versioned publication
  alongside the physics PDF from the tested source archive.
- A practical physics-manual chapter, shared reference-junction input, and a
  DMRG bond-refinement tutorial with ED comparisons for the same bath and checks of continuation from saved states.

### 0.3.0 (26 October 2026)

- Local impurity-orbital density matrices for QP and DMRG eigenstates, including
  incoherent mixtures, occupation probabilities, pairing coherences, and purity.
- Reusable QP workspaces and spectrum-coverage support.
- Cooperative interruption and cancellation of native QP operations.
- Additional published NRG comparison studies, reviewed scientific outputs, and
  validation coverage.

### 0.2.0 (2 October 2026)

The first public release.

- Shared finite Hamiltonians and physical conventions for QP and DMRG backends.
- QP calculations of several eigenstates of complex Hamiltonians, with orthonormal
  Rayleigh–Ritz refinement.
- `save_result`/`load_result` with one file format for both backends; earlier QP results
  remain readable.
- Companion wavefunction files that are never overwritten, and JSON results that replace
  the previous file only after a complete save.
- Strict validation of versioned input files and explicit convergence diagnostics.
- Independent finite-model tests, randomized (property-based) tests of the C++ extension,
  and CI checks for every pull request and every week.

See [interfaces](interfaces.md) for defaults and record-version compatibility,
and [numerics](numerics.md) for the distinct numerical convergence requirements.
