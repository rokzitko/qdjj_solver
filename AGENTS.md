# Repository guidance for coding agents

These instructions apply throughout this repository. Consult
[CONTRIBUTING.md](CONTRIBUTING.md) and the relevant files in `docs/` before
changing numerical behavior or public interfaces.

## Project and layout

`qdjj_solver` solves superconducting quantum-impurity models using a truncated
quasiparticle (QP) expansion or density-matrix renormalization group (DMRG).
Both backends consume the same physical Hamiltonian and conventions.

| Location | Responsibility |
|---|---|
| `src/qdjj_solver/common/` | Fermionic algebra, baths, models, sectors, configuration and records |
| `src/qdjj_solver/qp_solver/` | Finite-Fock-space QP diagonalization and projected operators |
| `cpp/core.cpp` | Native C++17/pybind11 basis and operator implementation |
| `src/qdjj_solver/dmrg_solver/` | TeNPy MPO/MPS calculations, excited states and checkpoints |
| `src/qdjj_solver/api.py`, `src/qdjj_solver/cli.py` | Public backend dispatch and unified command line |
| `src/qpsolver/` | Compatibility imports and legacy CLI; preserve this facade |
| `tests/` | Numerical, interface, persistence and packaging regressions |
| `examples/` | Small interface and workflow examples |
| `first_calculation/` | README calculation inputs, runner and reference outputs |
| `applications/` | Literature studies with inputs, references, outputs and reports |
| `docs/`, `tools/`, `.github/workflows/` | Documentation, validation utilities and CI |

Keep shared physics and serialization in `common`, rather than duplicating them
in the backends. Preserve lazy imports: the base/common API and QP functionality
must work without TeNPy or HDF5 dependencies. Follow the surrounding code's style
and maintain Python 3.11 compatibility.

## Development setup

Run commands from the repository root:

```sh
python -m pip install -e '.[dmrg,dev]'
```

Install the `plots` extra when generating figures. Reference-extraction tools may
have additional optional dependencies documented in their case README. Reinstall
the package after editing C++ or build configuration so Python uses the rebuilt
native extension. Runtime dependency lower bounds are in `pyproject.toml` and
`tools/minimum-runtime.txt`; keep them consistent when changing dependencies.

## Validation

Start with tests relevant to the change. For numerical or shared-API changes,
run the regular suite and the coverage gate:

```sh
python -B -m pytest --cov --cov-report=term-missing --cov-report=json
python tools/check_coverage.py reports/python-coverage.json
python -m ruff check src tests tools examples applications first_calculation
```

- Regular pytest runs exclude `slow`. Run applicable extended checks explicitly
  with `python -B -m pytest -m slow`; DMRG tests require the optional backend.
- CI requires **95% Python line coverage**. Branch/native coverage is diagnostic.
- For reproducible numerical runs, set `OPENBLAS_NUM_THREADS=1`,
  `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and, on macOS,
  `VECLIB_MAXIMUM_THREADS=1`.
- Use fixed seeds, explicit tolerances and independent physical checks: electron-
  space ED, analytic limits, symmetry identities, finite differences, or another
  backend. Explain numerical tolerances and convergence evidence.
- Documentation-only changes need relevant formatting/link and packaging checks;
  avoid rerunning large physics scans for prose edits.

For distribution/build changes, or changes to shipped assets:

```sh
python -m build
python tools/check_distribution.py dist/*.tar.gz dist/*.whl
python tools/check_wheel.py dist/*.whl --dmrg
```

The distribution checker has explicit source-root/file allowlists. Account for
new top-level files and required application assets. Applications ship with the
source distribution; wheels contain the installed solver packages. See
[docs/validation.md](docs/validation.md) for minimum-dependency and platform checks.

## Numerical conventions and pitfalls

- Preserve the definitions in [docs/models.md](docs/models.md). `bandwidth` is
  the normal-state **half-bandwidth**; bath weights represent the normal-state
  measure. Do not silently normalize coarse quadrature or surrogate weights.
- `reference_model(gamma=...)` uses **total** hybridization. Papers may instead
  quote per-lead couplings. `detuning=epsilon_d+u/2`; half filling is detuning zero.
- Phases are in radians. The reference junction uses `[-phi/2,+phi/2]` and
  `phi=phi_R-phi_L`. Positive local `field` raises the spin-up level.
- Python `solve(..., cutoff=None)` retains the full declared finite bath; the
  unified QP CLI defaults to cutoff two. Specify numerical settings explicitly
  in examples and research inputs.
- `compress=True` freezes discarded dark modes in their vacuum. Use
  `compress=False` for a full finite-bath spectrum and arbitrary electron-operator
  transformations. DMRG cannot recover modes removed by the model constructor.
- Find the global ground state by comparing relevant sectors. At finite field,
  consider both doublet spin projections. Transform physical electron operators
  into the Hamiltonian's canonical coordinates before evaluating transitions.
- The current is `(2e/hbar) * dE/dphi`. State its normalization explicitly;
  `e*Delta/hbar` and `2e*Delta/hbar` differ by two. Preserve energy-zero conventions.
- Distinguish eigensolver residual, QP cutoff/MPS bond convergence, and bath
  convergence. A small residual does not establish continuum accuracy.
- At exact degeneracy, verify the complete relevant subspace using dense or
  symmetry-resolved calculations. Individual eigenvectors and their expectations
  can be basis dependent; single-vector Lanczos may miss a repeated eigenvalue.

## Research applications and generated files

- Follow `applications/<case>/`: scientific `README.md`, `input/`, `reference/`,
  `output/`, and small runner/plot/extraction scripts. Reuse existing helpers.
- Identify the paper/DOI and exact figure or result. Document parameter mappings,
  units, reference provenance, extraction precision and the reproduction's scope.
- Retain compact reviewed numerical outputs, plots, explicit bath records and
  solver/environment provenance. Quantify literature discrepancies honestly.
- Keep routine regressions small and offline. Run full `paper`/`convergence`
  profiles explicitly; they may take minutes or hours. Respect solver resource
  limits and archive convergence evidence for the observables actually reported.
- Fresh runs and checkpoints belong in ignored `applications/<case>/runs/` or
  root `results/`/`runs/`. Generated test/build artifacts belong in ignored
  `reports/`, `build/` and `dist/`.
- **`APPLICATIONS.md` is a local research notebook. Keep it ignored and excluded
  from commits and distributions.** Essential scientific information must also
  appear in the public case reports.
- Update the application index, analyzer, bounded tests, CI selection and
  distribution checks when adding a case. Preserve user work and commit/push
  only when requested.

## Markdown and LaTeX: render correctly on GitHub

**Create and edit all `.md` files using GitHub-supported math syntax**, including
READMEs, application reports and other documentation.

- Use `$...$` for inline equations, for example `$E_D-E_S$` or `$\Gamma_L+\Gamma_R$`.
- For inline equations containing escaped punctuation such as `\_`, `\{`, `\}`
  or `\,`, use GitHub's protected dollar-and-backtick form: ``$`...`$``.
  For example, write ``$`\texttt{one\_body}`$`` or ``$`\{A,B\}`$``. Ordinary
  `$...$` spans can lose these backslashes during Markdown processing, causing
  MathJax errors, missing braces, or visible commas in place of thin spaces.
- Keep math separated from adjacent words: use `$k$-th`, not `$k$th`. Write
  ranges as `from $a$ to $b$` or as one math expression; GitHub can leave an
  immediately adjacent span such as `$a$–$b$` partly unrendered.
- For display equations, use a fenced **`math`** block (as in `docs/models.md`)
  or a `$$` block with each delimiter on its own line. Surround display blocks
  with blank lines.
- Avoid `\(...\)` and `\[...\]` delimiters in Markdown: they do not reliably
  render as mathematics on GitHub. Convert these in documentation sections you edit.
- Ordinary code spans and `tex`, `latex` or unlabeled code fences display literal
  source. Use them only when showing syntax or code, not for rendered equations.
- Use standard MathJax-compatible LaTeX, without document preambles, package
  imports or undeclared macros. Use `\mathrm{}` or `\text{}` for words in math.
  GitHub rejects `\operatorname`; use `\mathrm{asinh}` for the upright function
  name in `\mathrm{asinh}(D/\Delta)`, for example.
- Put multiline equations outside Markdown tables. In table cells, use inline
  math and `\vert`/`\lvert`/`\rvert` rather than literal pipe characters that
  interfere with table separators.
- Check delimiters, escaping, surrounding whitespace and links; preview on GitHub
  when available rather than assuming a local Markdown renderer is equivalent.
  Check that the intended symbols and spacing survive Markdown processing, as
  well as checking for MathJax errors.

Example of a GitHub-renderable display equation:

```math
I(\phi)=\frac{2e}{\hbar}\frac{\partial E_g(\phi)}{\partial\phi},
\qquad \Gamma=\Gamma_L+\Gamma_R.
```
