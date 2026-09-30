# Release artifacts and procedure

This is a **maintainer/developer guide**, not the current user installation
path. Users currently build from source following the
[root README](https://github.com/rokzitko/qdjj_solver/blob/main/README.md#install).
The PyPI publication steps below are for a future release; they do not imply
current PyPI availability.

## Versioned artifacts

`.github/workflows/release.yml` runs on `v*` tags and manual dispatch. A manual
run rehearses the full process and retains artifacts without publishing a
GitHub release, including when dispatched against a tag. Only a tag **push**
attaches the tested sdist and wheels to a GitHub release after all checks succeed.

The workflow:

1. Runs the shared CI suite, including minimum runtime dependencies, sanitizers,
   stress tests, and native coverage.
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

## Preparing a release

Before making the repository public, perform a **clean-history publication step**
from the reviewed source tree. The private development history includes local
machine paths in earlier research archives; deleting those files does not remove
them from older commits. Keep the development history privately, and validate
the clean public candidate and its freshly built distributions before publishing.
This is a pending publication step, not an instruction to rewrite history during
routine development.

Update the three version declarations together, refresh the copyable citation
in [references](references.md#software-citation), then check the declarations locally:

```sh
python -B tools/check_release.py --tag v0.2.0
```

Run a manual release-artifact workflow on the candidate commit and inspect all
job results. Confirm the public installation instructions and runnable examples
match the candidate. Create the matching tag on the verified commit when ready.
Before public release, verify that repository, documentation, and artifact links
are accessible without authentication.

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

## First public release notes

- Shared finite Hamiltonians and physical conventions for QP and DMRG backends.
- Complex multiroot QP calculations with orthonormal Rayleigh–Ritz refinement.
- Unified common-format `save_result`/`load_result`, with legacy QP compatibility.
- Immutable companion state artifacts and failure-atomic JSON publication.
- Strict versioned configuration validation and explicit convergence diagnostics.
- Independent finite-model tests, property-based native checks, and tiered CI.

See [interfaces](interfaces.md) for defaults and record-version compatibility,
and [numerics](numerics.md) for the distinct numerical convergence requirements.
