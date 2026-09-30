import json
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import (Hamiltonian, Sector, cosh_grid, create, number,
                         reference_model, solve)
from qdjj_solver.common.io import fingerprint, model_from_config, source_provenance, write_json


def test_shared_types_and_legacy_results(tmp_path):
    import qpsolver
    from qpsolver.io import calculation_from_config, save_result
    from qdjj_solver.common import Hamiltonian as SharedHamiltonian
    from qdjj_solver.qp_solver import Hamiltonian as QPHamiltonian
    assert SharedHamiltonian is QPHamiltonian is qpsolver.Hamiltonian
    import qpsolver.models
    import qpsolver._core
    from qdjj_solver.qp_solver import _core
    assert qpsolver.models.Hamiltonian is SharedHamiltonian
    assert qpsolver._core is _core
    config = dict(model=dict(geometry="single"), bath=dict(pairs=1, bandwidth=5), cutoff=2)
    result = qpsolver.solve(*calculation_from_config(config))
    path = tmp_path/"legacy.json"
    save_result(result, path, save_vectors=True)
    record = json.loads(path.read_text())
    assert record["format"] == "bcs-qp-eigenstates"
    assert record["format_version"] == 1
    restored = qpsolver.solve(Hamiltonian.from_record(record["hamiltonian"]), 2)
    assert_allclose(restored.energies, result.energies, atol=1e-13)
    assert restored.hamiltonian.fingerprint() == result.hamiltonian.fingerprint()


def test_common_import_does_not_load_numerical_backends():
    command = [sys.executable, "-B", "-c", "import qdjj_solver, sys; "
               "assert 'tenpy' not in sys.modules; assert 'qdjj_solver.qp_solver._core' not in sys.modules; "
               "import qpsolver; assert 'tenpy' not in sys.modules"]
    subprocess.run(command, check=True, capture_output=True, text=True)


def test_physical_maps_and_coordinate_identity():
    bath = cosh_grid(pairs=1, bandwidth=4)
    h = reference_model(bath, compress=False)
    h2 = reference_model(bath, u=3.1, compress=False)
    assert h.fingerprint() != h2.fingerprint()
    assert h.coordinate_fingerprint() == h2.coordinate_fingerprint()
    charge_difference = h.physical_operator(number(0)+number(1))-h.observables["impurity_charge"]
    assert max(map(abs, charge_difference.terms.values()), default=0.) < 1e-14
    physical_spin = (number(0)-number(1))/2
    difference = h.physical_operator(physical_spin)-h.observables["impurity_spin_z"]
    assert max(map(abs, difference.terms.values()), default=0.) < 1e-14
    h_compressed = reference_model(bath, compress=True)
    assert h.coordinate_fingerprint() != h_compressed.coordinate_fingerprint()
    with pytest.raises(ValueError, match="uncompressed"):
        h_compressed.physical_operator(physical_spin)


def test_atomic_json_and_provenance(tmp_path):
    path = tmp_path/"state.json"
    write_json(path, {"value": np.array([1+2j])})
    previous = path.read_bytes()
    with pytest.raises(ValueError):
        write_json(path, {"value": float("nan")})
    assert path.read_bytes() == previous
    assert fingerprint({"b": np.int64(3), "a": 2}) == fingerprint({"a": 2, "b": 3})
    provenance = source_provenance()
    assert "common/problem.py" in provenance["files"]
    assert "compat/qpsolver/__init__.py" in provenance["files"]
    assert all(not name.startswith("/") for name in provenance["files"])


def test_unified_cli_and_legacy_cli_outside_checkout(tmp_path):
    config = tmp_path/"model.json"
    write_json(config, dict(model=dict(geometry="single"), bath=dict(pairs=1), cutoff=2))
    for module, expected in (("qdjj_solver.cli", "qdjj-eigenstates"), ("qpsolver.cli", "bcs-qp-eigenstates")):
        output = tmp_path/f"{module}.json"
        completed = subprocess.run([sys.executable, "-B", "-m", module, "solve", str(config),
                                    "--output", str(output)], cwd=tmp_path,
                                   text=True, capture_output=True, check=True)
        assert "State 0:" in completed.stdout
        assert json.loads(output.read_text())["format"] == expected


def test_dmrg_checkpoint_and_continuation(tmp_path):
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    from qdjj_solver.dmrg_solver import SolverOptions
    from qdjj_solver.dmrg_solver.convergence import compare_states
    from qdjj_solver.dmrg_solver.io import load_result, save_result
    h = reference_model(cosh_grid(pairs=1, bandwidth=3), compress=False)
    options = SolverOptions(eigenpairs=2, chi_max=16, require_convergence=True,
                            observables=("impurity_spin_z", "impurity_charge"))
    result = solve(h, backend="dmrg", options=options, sector=Sector(1, 1, 1))
    exact = solve(h, backend="qp", sector=Sector(1, 1, 1), options={"eigenpairs": 2})
    assert_allclose(result.energies, exact.energies, atol=1e-10)
    path = tmp_path/"state.json"
    save_result(result, path, save_states=True)
    restored = load_result(path)
    assert_allclose(restored.energies, result.energies, atol=1e-14)
    assert_allclose(abs(restored.overlaps(result)), np.eye(2), atol=1e-10)
    assert_allclose(restored.matrix_elements(h.observables["impurity_spin_z"]),
                    result.matrix_elements(h.observables["impurity_spin_z"]), atol=1e-12)
    continued = solve(h, backend="dmrg", options=options, sector=Sector(1, 1, 1), initial=restored)
    comparison = compare_states(result, continued)
    assert_allclose(comparison["energy_changes"], 0., atol=1e-11)
    record = json.loads(path.read_text())
    record["energies"][0] += 1.
    write_json(path, record)
    with pytest.raises(ValueError, match="do not match"):
        load_result(path)


def test_unified_backend_selection_preserves_the_finite_model(tmp_path):
    pytest.importorskip("tenpy")
    records = []
    for backend in ("qp", "dmrg"):
        config, output = tmp_path/f"{backend}-input.json", tmp_path/f"{backend}-result.json"
        write_json(config, dict(backend=backend, cutoff=None,
                                model=dict(u=1.6, gamma=.3, phi=.8),
                                bath=dict(pairs=1, bandwidth=3.)))
        subprocess.run([sys.executable, "-B", "-m", "qdjj_solver.cli", "solve", str(config),
                        "--output", str(output)], cwd=tmp_path, check=True, capture_output=True, text=True)
        records.append(json.loads(output.read_text()))
    assert records[0]["hamiltonian_sha256"] == records[1]["hamiltonian_sha256"]
    assert not records[0]["hamiltonian"]["metadata"]["paired_mode_compression"]
    assert_allclose(records[0]["energies"], records[1]["energies"], atol=2e-10)


def test_complex_config_and_backend_options():
    config = dict(model=dict(kind="multi-orbital",
                            impurity=dict(orbitals=1, terms=(number(0)+number(1)).to_records()),
                            reservoirs=[dict(kind="discrete", xi=[0.], weights=[1.], delta=1., bandwidth=2.)],
                            tunneling=[{"real": [[.1, 0.], [0., .1]], "imag": [[0., .02], [.02, 0.]]}]))
    h = model_from_config(config)
    assert h.operator.hermiticity_error() < 1e-14
    with pytest.raises(ValueError, match="backend"):
        solve(h, backend="unknown")
    with pytest.raises(ValueError, match="parity"):
        Hamiltonian(create(0)+create(0).dagger(), 2, (1, -1))
