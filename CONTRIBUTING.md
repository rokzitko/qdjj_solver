# Contributing

First satisfy the [source-build prerequisites](https://github.com/rokzitko/qdjj_solver/blob/main/README.md#install),
then install from the repository root with `python -m pip install -e '.[dmrg,dev]'`.
Keep physical models, operator algebra, coordinate conventions, and codecs in
`common`. Backends consume that contract rather than maintaining separate
physical conventions. Preserve the installed `qpsolver` façade.

For numerical changes, use independent comparisons appropriate to the affected
physics: tensor-product electron ED, native projected actions, analytic limits,
symmetry identities, or physical finite differences. Include convergence
evidence rather than merely copying implementation outputs into fixtures.
A residual is not a cutoff or continuum error estimate.

```sh
python -B -m pytest --cov --cov-report=term-missing --cov-report=json
python tools/check_coverage.py reports/python-coverage.json
python -m ruff check src tests tools examples applications first_calculation large_Gamma NRG_comparisons
python -m build
python tools/check_distribution.py dist/*.tar.gz dist/*.whl
python tools/check_wheel.py dist/*.whl --dmrg
```

The wheel checker accepts **one wheel**. If `dist/` contains multiple wheels,
replace `dist/*.whl` in the last command with the path to one selected wheel
compatible with your Python interpreter and platform. Keep `--dmrg` to check
both backends and HDF5 persistence.

The regular suite excludes `slow` tests. Run those explicitly with
`python -B -m pytest -m slow`. CI requires 95% Python line coverage; branch and
native coverage are diagnostics. Prefer independent numerical oracles and
focused regression tests over assertions about private implementation details.
See [validation](docs/validation.md) for the CI jobs and native instrumentation.
Keep `tools/minimum-runtime.txt` synchronized with supported runtime lower bounds.
Tagged builds and version agreement are described in [releases](docs/releases.md).

Keep interface examples small. Literature reproductions belong in `applications/`,
with explicit inputs, compact reference/output files, source-data provenance,
and a scientific README. Add bounded regressions under `tests/`; full paper scans
are explicit application commands. Record bath and many-body convergence
separately, and keep fresh runs/checkpoints in ignored `runs/` directories.
