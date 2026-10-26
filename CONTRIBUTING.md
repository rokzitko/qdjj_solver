# Contributing

First satisfy the [source-build prerequisites](README.md#install), then install from the
repository root with `python -m pip install -e '.[dmrg,dev,docs]'`. The
[architecture and extension guide](docs/developer.md) describes the layout, the design rules
and how to add models, observables, probes, drives, reservoirs and options. In short: keep
physical models, operator algebra, coordinate conventions, and the reading and writing of
records in `common`; the backends use these shared definitions and keep no physical
conventions of their own; preserve the installed `qpsolver` façade.

For numerical changes, use independent comparisons appropriate to the affected
physics: tensor-product electron ED, the projected operator actions of the C++
extension, analytic limits, symmetry identities, or physical finite differences.
Include convergence evidence; do not merely copy the code's own output into the
reference values of a test. A residual is not an estimate of the cutoff or
continuum error.

Before submitting, run the [local checks](docs/validation.md#local-checks): the test suite
with the 95% line-coverage gate, `ruff`, the distribution checks and, when documented
behaviour or documentation changes, the documentation checks. The regular suite excludes
the `slow`, `docs` and `tutorial` tests; run them explicitly with `-m slow`, `-m docs` or
`-m tutorial`. Keep `tools/minimum-runtime.txt` synchronized with the supported runtime lower
bounds in `pyproject.toml`. Tagged builds and version agreement are described in
[releases](docs/releases.md).

**Documentation is part of every change.** A change of features or behaviour updates the
physics manual, the user manual and a tutorial, or the pull request states why one of them
is unaffected. Reference tables are generated from the code with `python tools/gen_docs.py`.
The [documentation policy](docs/documentation.md) explains which topics each document covers and how
the checks keep them in sync.

Keep interface examples small. Literature reproductions belong in `applications/`,
with explicit inputs, compact reference/output files, source-data provenance,
and a scientific README. Add small regression tests under `tests/`; full literature-reproduction scans
are explicit application commands. Record bath and many-body convergence
separately, and keep fresh runs and checkpoints in ignored `runs/` directories.
