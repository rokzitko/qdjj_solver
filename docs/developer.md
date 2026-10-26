# Architecture and extension guide

This page is for developers who change or extend the package. The physical definitions
are defined in the physics manual and summarized in the
[conventions reference](reference/conventions.md); the test and CI policy is in
[validation](validation.md), and the documentation rules are in
[documentation](documentation.md).

## Layout

| Location | Responsibility |
|---|---|
| `src/qdjj_solver/common/` | Backend-independent physics: fermionic algebra (`algebra.py`), reservoirs (`baths.py`), models and observables (`models.py`), Hamiltonians and sectors (`problem.py`), periodic drives (`floquet.py`), probes (`probes.py`), spectra (`spectra.py`), configuration files and records (`io.py`) |
| `src/qdjj_solver/qp_solver/` | Finite-Fock-space QP diagonalization (`solver.py`, `basis.py`), Floquet methods (`floquet.py`), Floquet matrix elements (`transitions.py`), persistence and the legacy command line |
| `cpp/core.cpp` | Native C++17/pybind11 Fock basis and operator action, used by `qp_solver` |
| `src/qdjj_solver/dmrg_solver/` | TeNPy MPO/MPS calculations, excited states, residuals and checkpoints |
| `src/qdjj_solver/api.py`, `cli.py` | Backend dispatch (`solve`, `save_result`, `load_result`) and the unified command line |
| `src/qpsolver/` | Compatibility imports and the legacy `qpsolver` command; keep this façade working |
| `tests/` | Numerical, interface, persistence, packaging and documentation checks |
| `examples/` | Tutorials (`examples/NN_topic/`) and small interface examples |
| `applications/`, `large_Gamma/`, `NRG_comparisons/` | Literature reproductions and research studies |
| `docs/`, `docs/manual/`, `tools/` | User and developer manual, physics manual, generators and checks |

## Design rules

- **Shared physical definitions.** Models, operator algebra, coordinate conventions and serialization routines
  belong in `common`. The backends use the same `Hamiltonian` and do not maintain separate
  physical conventions.
- **Lazy optional dependencies.** The base package and the QP backend must work without
  TeNPy or HDF5; import them only inside `dmrg_solver`, and dispatch to the DMRG backend
  lazily in `api.py`.
- **Canonical coordinates.** Every operator, observable and probe acts in the canonical
  (Bogoliubov) coordinates of its model. Physical electron operators are transformed with
  `Hamiltonian.physical_operator`, which requires an uncompressed model; results with
  different `coordinate_fingerprint()` values cannot be combined.
- **Explicit numerical settings.** Options are frozen dataclasses with validated fields and
  a docstring per field; the [option reference](reference/options.md) is generated from them.
- **Python 3.11 compatibility** and the style of the surrounding code; `ruff` runs in CI.

## Numerical conventions and pitfalls

- `bandwidth` is the normal-state **half**-bandwidth; bath weights represent the
  normal-state measure. Do not silently normalize coarse quadrature or surrogate weights.
- `reference_model(gamma=...)` uses the **total** hybridization; papers may quote per-lead
  couplings. `detuning=epsilon_d+u/2`, so half filling is detuning zero.
- Phases are in radians; the reference junction uses $[-\phi/2,+\phi/2]$ with
  $\phi=\phi_R-\phi_L$. A positive local `field` raises the spin-up level.
- Python `solve(..., cutoff=None)` keeps the full declared finite bath, while the
  `qdjj-solver` command defaults to cutoff two.
- `compress=True` freezes discarded dark modes in their vacuum; use `compress=False` for a
  full finite-bath spectrum, for electron operators and for every drive.
- Find the global ground state by comparing the relevant sectors; at finite field, include
  both doublet spin projections.
- The current is $`(2e/\hbar)\,\partial E/\partial\phi`$, positive from the left to the right
  lead; state its normalization explicitly. The bias voltage is $V=V_L-V_R$, so
  $\dot\phi=2eV/\hbar$.
- Distinguish the eigensolver residual, QP-cutoff or bond-dimension convergence, and bath
  convergence: a small residual does not establish continuum accuracy.
- At exact degeneracy compare complete multiplets; individual eigenvectors are basis
  dependent and single-vector Lanczos may miss a repeated eigenvalue.

## Extending the package

Each extension needs an independent test ([validation](validation.md)) and documentation
in all three places ([documentation](documentation.md)). Typical cases:

- **A model or impurity.** Build it with `Impurity` (or `Impurity.from_integrals`) and
  `make_model`; declare observables in canonical coordinates and keep the
  electron-to-canonical map in `metadata["coordinates"]`. Test against the independent
  tensor-product electron Hamiltonian of `tests/electron_oracle.py`.
- **A built-in observable.** Add it to the model and describe it in `OBSERVABLE_DOCS`
  (`common/models.py`); `tests/test_docs.py` requires the descriptions to cover exactly the
  built-in observables.
- **A probe.** In `common/probes.py`, add the name to `PROBES` and its description to
  `PROBE_DOCS`, in the same order. Map an observable probe to its model observable in
  `OBSERVABLE_PROBES`; for an electron operator, add the name to `ELECTRON_PROBES` and the
  operator to `electron_probe`.
- **A drivable parameter.** Any builder keyword can be driven with `ParameterDrive`, provided
  the canonical coordinates stay fixed during the period; `PeriodicHamiltonian` rejects
  drives that move them.
- **A reservoir representation.** Return a `DiscreteBath` whose weights represent the
  normal-state measure, and document its parameters and convergence in
  [bath representations](bath_representations.md).
- **An option or result field.** Add the dataclass field with a docstring (plain text and
  double-backtick literals), validate it in `__post_init__`, then run
  `python tools/gen_docs.py` to regenerate the reference pages.
- **A user-facing error message.** Add an entry to `docs/troubleshooting.toml`;
  the generator checks that the quoted message occurs in the code.

## Native extension

`cpp/core.cpp` enumerates Fock bases and applies operators. Reinstall the package
(`python -m pip install -e '.[dmrg,dev]'`) after changing it or the build configuration, so
that Python uses the rebuilt extension. The sanitizer and native-coverage CI jobs are
described in [validation](validation.md#native-instrumentation).
