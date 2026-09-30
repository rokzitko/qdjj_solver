"""Portable records and failure-atomic publication, including real state artifacts."""

from copy import deepcopy
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import pickle
from unittest.mock import Mock

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from qdjj_solver import Hamiltonian, Sector, create, load_result, number, save_result, solve
from qdjj_solver.common import baths as bath_module, io
from qdjj_solver.qp_solver import io as qp_io


def _paired_hamiltonian():
    pair = (.2+.3j)*create(0)*create(1)
    return Hamiltonian(.7*number(0)+1.1*number(1)+pair+pair.dagger(), 2, (1, -1),
                       observables={"charge": number(0)+number(1), "pair": create(0)*create(1)},
                       metadata={"samples": np.array([1., 2.]), "phase": np.array([1+2j])})


@pytest.fixture(scope="module")
def qp_result():
    return solve(_paired_hamiltonian(), sector=Sector(0, 0), options={"eigenpairs": 2})


@pytest.fixture(scope="module")
def mps_result():
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    local = _paired_hamiltonian()
    h = Hamiltonian(local.operator+sum((4+i)*number(i) for i in range(2, 6)), 2, (1, -1)*3,
                    observables=local.observables, metadata=local.metadata)
    return solve(h, sector=Sector(0, 0), backend="dmrg",
                 options=dict(eigenpairs=2, group_size=2, mode_order=(1, 0, 2, 3, 4, 5), seed_trials=1,
                              chi_max=8, require_convergence=True))


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_numpy_json_encoding_and_complex_decoding():
    source = {7: (np.int64(3), np.float32(1.25), np.bool_(True), None),
              "real": np.arange(6).reshape(2, 3), "complex": np.array([[1+2j, 3-4j]]),
              "scalar": np.complex128(2-3j), "plain": 4+5j}
    encoded = io.json_value(source)
    assert json.loads(json.dumps(encoded))["7"] == [3, 1.25, True, None]
    assert_array_equal(encoded["real"], source["real"])
    assert_allclose(io.complex_array(encoded["complex"]), source["complex"], rtol=0, atol=0)
    assert io.complex_array(encoded["scalar"]) == 2-3j
    assert encoded["plain"] == {"real": 4., "imag": 5.}
    assert_array_equal(io.complex_array([1, 2]), [1+0j, 2+0j])
    assert io.complex_array({"real": [], "imag": []}).shape == (0,)
    assert source["real"].shape == (2, 3)


@pytest.mark.parametrize("value, message", [
    ({"real": [1]}, "real and imag fields"),
    ({"real": [1, 2], "imag": [0]}, "shapes differ"),
])
def test_malformed_complex_array(value, message):
    with pytest.raises(ValueError, match=message):
        io.complex_array(value)


def test_fingerprints_are_canonical_and_reject_nonfinite():
    expected = hashlib.sha256(b'{"a":[1,2],"b":3}').hexdigest()
    assert io.fingerprint({"b": np.int64(3), "a": np.array([1, 2])}) == expected
    assert io.fingerprint({"a": (1, 2), "b": 3}) == expected
    assert io.fingerprint({"a": [2, 1], "b": 3}) != expected
    with pytest.raises(ValueError, match="Out of range"):
        io.fingerprint({"value": np.nan})


def test_atomic_json_failure_preserves_previous_record(tmp_path, monkeypatch):
    path = tmp_path/"nested"/"record.json"
    io.write_json(path, {"old": np.arange(2)})
    before = path.read_bytes()

    def fail(*args, **kwargs):
        raise OSError("injected publication failure")

    monkeypatch.setattr(io.os, "replace", fail)
    with pytest.raises(OSError, match="injected publication failure"):
        io.write_json(path, {"new": 3})
    assert path.read_bytes() == before
    assert list(path.parent.iterdir()) == [path]


def test_atomic_json_replacement_and_portable_text(tmp_path, monkeypatch):
    path = tmp_path/"record.json"
    io.write_json(path, {"old": True})
    handles = []
    temporary_file, replace_file = io.tempfile.NamedTemporaryFile, io.os.replace

    def track_temporary(*args, **kwargs):
        handle = temporary_file(*args, **kwargs)
        handles.append(handle)
        return handle

    def replace_closed(*args, **kwargs):
        assert handles and all(handle.closed for handle in handles)
        return replace_file(*args, **kwargs)

    monkeypatch.setattr(io.tempfile, "NamedTemporaryFile", track_temporary)
    monkeypatch.setattr(io.os, "replace", replace_closed)
    io.write_json(path, {"z": np.complex128(1+2j), "a": np.array([3., 4.])})
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    assert text.index('"a"') < text.index('"z"')
    assert json.loads(text) == {"a": [3., 4.], "z": {"real": 1., "imag": 2.}}
    assert list(tmp_path.iterdir()) == [path]


def test_numerical_environment_handles_missing_packages(monkeypatch):
    def version(name):
        if name in ("physics-tenpy", "h5py", "scipy"):
            raise io.metadata.PackageNotFoundError(name)
        return f"test-{name}"

    monkeypatch.setattr(io.metadata, "version", version)
    environment = io.numerical_environment()
    assert environment["qdjj_solver"] == environment["qpsolver"]
    assert environment["numpy"] == np.__version__
    assert environment["python"]
    assert environment["scipy"] is None
    assert "h5py" not in environment["packages"]


def test_bath_and_model_config_routes_preserve_input():
    default = io.bath_from_config({"pairs": 1, "bandwidth": 3.})
    assert default.levels == 2
    assert default.metadata["kind"] == "cosh-grid"
    fitted = io.bath_from_config(dict(kind="surrogate", levels=1, bandwidth=3.,
                                      frequency_points=12, starts=1))
    assert fitted.levels == 1
    assert fitted.metadata["kind"] == "surrogate"
    assert fitted.metadata["max_nfev"] == 5000
    discrete = io.bath_from_config(dict(kind="discrete", xi=[-.5, .2], weights=[.3, .7]))
    assert_array_equal(discrete.xi, [-.5, .2])
    with pytest.raises(ValueError, match="unknown bath"):
        io.bath_from_config({"kind": "unrecognized"})
    config = {"bath": {"pairs": 1, "bandwidth": 3.}}
    original = deepcopy(config)
    compressed = io.model_from_config(config)
    full = io.model_from_config(config, default_compress=False)
    assert config == original
    assert compressed.metadata["paired_mode_compression"]
    assert not full.metadata["paired_mode_compression"]
    explicit = io.model_from_config(dict(config, model={"compress": False}))
    assert explicit.fingerprint() == full.fingerprint()
    with pytest.raises(ValueError, match="model kind"):
        io.model_from_config({"model": {"kind": "other"}})


def test_surrogate_config_forwards_budget_and_roundtrips_diagnostics(tmp_path, monkeypatch):
    config = dict(kind="surrogate", levels=1, bandwidth=3., frequency_points=12,
                  starts=1, max_nfev=73)
    original = deepcopy(config)
    optimizer = Mock(wraps=bath_module.least_squares)
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    bath = io.bath_from_config(config)
    assert config == original
    optimizer.assert_called_once()
    assert optimizer.call_args.kwargs["max_nfev"] == 73
    metadata = bath.metadata
    assert metadata["max_nfev"] == 73
    assert metadata["selected_start"] == 0
    trial, = metadata["optimizer_trials"]
    assert trial["success"] is True and trial["status"] > 0
    assert isinstance(trial["message"], str) and trial["message"]
    assert 0 < trial["evaluations"] == metadata["evaluations"] <= 73
    assert np.isfinite(trial["optimality"]) and trial["optimality"] == metadata["fit_optimality"]
    assert np.isfinite(trial["cost"]) and trial["cost"] == metadata["cost"]
    record = dict(kind="discrete", **bath.record())
    encoded = json.loads(json.dumps(record, allow_nan=False))
    path = tmp_path/"surrogate.json"
    io.write_json(path, record)
    assert json.loads(path.read_text()) == encoded
    restored = io.bath_from_config(encoded)
    assert restored.record() == bath.record()


@pytest.mark.parametrize("max_nfev", [0, True, 1.5])
def test_surrogate_config_rejects_invalid_budget_with_path(max_nfev):
    with pytest.raises(ValueError, match="bath: max_nfev must be a positive integer"):
        io.bath_from_config(dict(kind="surrogate", levels=1, frequency_points=10,
                                 starts=1, max_nfev=max_nfev))


def test_multiorbital_complex_config_matches_physical_electron_matrix():
    from electron_oracle import electron_hamiltonian, electron_operators, evaluate
    from qdjj_solver.common.baths import DiscreteBath
    from qdjj_solver.common.models import Impurity

    impurity = Impurity(1, .4*number(0)+.6*number(1), {"charge": number(0)+number(1)}, {"tag": "test"})
    baths = [DiscreteBath([-.2], [.3], bandwidth=2.), DiscreteBath([.4], [.7], bandwidth=2.)]
    tunneling = [np.array([[.1, .02j], [.03j, .2]]), np.diag([.13, .17])]
    direct = np.array([[.05j, .01], [-.02j, -.05j]])
    config = io.json_value(dict(model=dict(
        kind="multi-orbital", impurity=dict(orbitals=1, terms=impurity.operator.to_records(),
                                           observables={"charge": impurity.observables["charge"].to_records()},
                                           metadata=impurity.metadata),
        reservoirs=[dict(kind="discrete", **bath.record()) for bath in baths],
        tunneling=tunneling, phases=[-.2, .3],
        direct=[dict(leads=[0, 1], matrix=direct)])))
    before = deepcopy(config)
    h = io.model_from_config(config)
    physical = electron_hamiltonian(impurity, baths, tunneling, [-.2, .3], {(0, 1): direct})
    represented = evaluate(h.operator, electron_operators(6))
    assert_allclose(np.linalg.eigvalsh(represented.toarray()), np.linalg.eigvalsh(physical.toarray()),
                    rtol=0, atol=2e-14)
    assert "charge" in h.observables
    assert config == before


@pytest.mark.parametrize("changes,message", [
    ({"impurity": []}, "model.impurity must be an object"),
    ({"impurity": {"orbitals": 1, "terms": [], "typo": 1}}, "unknown model.impurity keys"),
    ({"impurity": {"orbitals": 1}}, "missing model.impurity keys"),
    ({"impurity": {"orbitals": 1, "terms": {}}}, "model.impurity.terms must be an array"),
    ({"impurity": {"orbitals": 1, "terms": [dict(operators=[1, -1], real=1, imag=0, typo=1)]}},
     r"unknown model.impurity.terms\[0\] keys"),
    ({"impurity": {"orbitals": 1, "terms": [], "observables": []}}, "model.impurity.observables must be an object"),
    ({"reservoirs": {}}, "model.reservoirs must be an array"),
    ({"reservoirs": [[]]}, r"model.reservoirs\[0\] must be an object"),
    ({"reservoirs": [{"pairs": 1, "bandwith": 3}]}, r"unknown model.reservoirs\[0\] keys"),
    ({"tunneling": {}}, "model.tunneling must be an array"),
    ({"tunneling": [{"real": [[0, 0], [0, 0]]}]}, r"model.tunneling\[0\].*real and imag"),
    ({"direct": {}}, "model.direct must be an array"),
    ({"direct": [{"leads": [0, 1]}]}, r"missing model.direct\[0\] keys"),
    ({"direct": [{"leads": [0., 1.], "matrix": [[0, 0], [0, 0]]}]}, r"model.direct\[0\].leads"),
    ({"direct": [{"leads": [0, 1], "matrix": [[0, 0], [0, 0]]}]*2}, "duplicates contact"),
])
def test_nested_model_configuration_errors_include_the_field_path(changes, message):
    model = dict(kind="multi-orbital", impurity=dict(orbitals=1, terms=[]),
                 reservoirs=[dict(pairs=1)], tunneling=[[[.1, 0], [0, .1]]])
    model.update(changes)
    with pytest.raises(ValueError, match=message):
        io.model_from_config({"model": model})


@pytest.mark.parametrize("observable", [create(2), create(2).dagger(), number(2), number(20)])
def test_undeclared_impurity_observable_fails_before_substitution(tmp_path, capsys, monkeypatch, observable):
    from qdjj_solver import Impurity, cli

    message = "impurity observable 'invalid' acts on an undeclared fermionic mode"
    with pytest.raises(ValueError, match=message):
        Impurity(1, number(0), observables={"invalid": observable})
    config = dict(model=dict(kind="multi-orbital",
                            impurity=dict(orbitals=1, terms=number(0).to_records(),
                                          observables={"invalid": observable.to_records()}),
                            reservoirs=[dict(pairs=1)], tunneling=[[[.1, 0], [0, .1]]]))
    with pytest.raises(ValueError, match=f"model.impurity: {message}"):
        io.model_from_config(config)
    solver = Mock(side_effect=AssertionError("invalid observable reached the solver"))
    monkeypatch.setattr(cli, "solve", solver)
    path, output = tmp_path/"model.json", tmp_path/"result.json"
    io.write_json(path, config)
    assert cli.main(["solve", str(path), "--output", str(output)]) == 2
    captured = capsys.readouterr()
    assert f"model.impurity: {message}" in captured.err and "Traceback" not in captured.err
    assert not output.exists()
    solver.assert_not_called()


def test_impurity_observable_validation_allows_odd_and_nonhermitian_operators():
    from qdjj_solver import Impurity

    observables = {"add_down": create(1), "pair": create(0)*create(1)}
    impurity = Impurity(1, number(0), observables=observables)
    for name, operator in observables.items():
        assert impurity.observables[name].terms == operator.terms


def test_explicit_hamiltonian_config_rejects_ignored_fields():
    record = dict(kind="hamiltonian", **_paired_hamiltonian().record(), misspelled_field=True)
    with pytest.raises(ValueError, match="unknown model keys.*misspelled_field"):
        io.model_from_config({"model": record})


def test_hamiltonian_config_and_ndarray_metadata_roundtrip(tmp_path):
    h = _paired_hamiltonian()
    path = tmp_path/"model.json"
    io.write_json(path, {"model": {"kind": "hamiltonian", **h.record()}})
    restored = io.model_from_config(json.loads(path.read_text()))
    assert restored.fingerprint() == h.fingerprint()
    assert_array_equal(restored.metadata["samples"], h.metadata["samples"])
    assert_array_equal(io.complex_array(restored.metadata["phase"]), h.metadata["phase"])
    assert restored.operator.terms == h.operator.terms


@pytest.mark.parametrize("config, message", [
    ({"unused": 1}, "unknown configuration keys.*unused"),
    ({"format_version": 2}, "unsupported configuration version"),
])
def test_qp_config_rejects_unknown_schema(config, message):
    with pytest.raises(ValueError, match=message):
        qp_io.calculation_from_config(config)


def test_qp_config_defaults_and_explicit_options():
    config = dict(model={"geometry": "single"}, bath={"pairs": 1}, cutoff=None,
                  sector={"parity": 0, "twice_sz": 0}, solver={"eigenpairs": 2, "method": "dense"})
    h, cutoff, sector, options = qp_io.calculation_from_config(config)
    assert cutoff is None and sector == Sector(0, 0)
    assert options.eigenpairs == 2 and options.method == "dense"
    assert h.metadata["geometry"] == "single-reservoir"
    _, cutoff, sector, options = qp_io.calculation_from_config({"bath": {"pairs": 1}})
    assert cutoff == 2 and sector == Sector() and options.eigenpairs == 1


def test_result_envelopes_and_missing_hamiltonian(qp_result):
    import qpsolver
    from qdjj_solver.qp_solver import _core

    for included in (False, True):
        legacy = qp_io.result_record(qp_result, include_hamiltonian=included)
        common = io.result_record(qp_result, backend="qp", include_hamiltonian=included)
        assert ("hamiltonian" in legacy) is included
        assert ("hamiltonian" in common) is included
        assert common["hamiltonian_sha256"] == qp_result.hamiltonian.fingerprint()
        assert common["backend"] == "qp"
        # Derive paths from imported modules, including for a wheel outside the checkout.
        files = common["implementation"]["files"]
        assert files["common/io.py"] == _sha(Path(io.__file__))
        assert files["compat/qpsolver/__init__.py"] == _sha(Path(qpsolver.__file__))
        assert files[f"qp_solver/{Path(_core.__file__).name}"] == _sha(Path(_core.__file__))
        assert all(not Path(name).is_absolute() for name in files)
        canonical = json.dumps(files, sort_keys=True, separators=(",", ":")).encode("utf-8")
        assert common["implementation"]["sha256"] == hashlib.sha256(canonical).hexdigest()
    detached = replace(qp_result, hamiltonian=None)
    assert "hamiltonian" not in qp_io.result_record(detached)
    with pytest.raises(ValueError, match="no Hamiltonian"):
        qp_io.result_record(detached, include_hamiltonian=True)


def test_qp_npz_roundtrip_hash_and_independent_eigen_equation(tmp_path, qp_result):
    path = tmp_path/"nested"/"result.json"
    qp_io.save_result(qp_result, path, save_vectors=True, common_format=True)
    record = json.loads(path.read_text())
    assert record["format"] == "qdjj-eigenstates"
    assert_array_equal(io.complex_array(record["hamiltonian"]["metadata"]["phase"]), [1+2j])
    artifact = record["state_amplitudes"]
    assert Path(artifact["file"]).name == artifact["file"]
    npz = path.parent/artifact["file"]
    assert artifact["sha256"] == _sha(npz)
    assert record["state_artifact"] == dict(kind="fock-amplitudes-npz", format_version=1, **artifact)
    with np.load(npz, allow_pickle=False) as arrays:
        assert set(arrays.files) == {"energies", "vectors", "occupations", "qp_counts", "result_sha256"}
        assert_array_equal(arrays["occupations"], [[0], [3]])
        assert_array_equal(arrays["qp_counts"], [0, 0])
        assert_array_equal(arrays["energies"], qp_result.energies)
        matrix = np.array([[0., .2-.3j], [.2+.3j, 1.8]])
        assert_allclose(matrix @ arrays["vectors"], arrays["vectors"]*arrays["energies"], rtol=0, atol=1e-14)
        assert_allclose(arrays["vectors"].conj().T @ arrays["vectors"], np.eye(2), rtol=0, atol=1e-14)


def test_qp_save_rejects_non_json_filename(tmp_path, qp_result):
    with pytest.raises(ValueError, match="end in .json"):
        qp_io.save_result(qp_result, tmp_path/"result.npz", True)
    assert not list(tmp_path.iterdir())


def test_qp_failed_overwrite_preserves_old_pair(tmp_path, qp_result, monkeypatch):
    path = tmp_path/"result.json"
    qp_io.save_result(qp_result, path, True, common_format=True)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    changed = replace(qp_result, energies=qp_result.energies+.1)

    def fail(*args, **kwargs):
        raise OSError("injected publication failure")

    # Fail after the replacement NPZ was written, at JSON publication.
    monkeypatch.setattr(qp_io, "write_json", fail)
    with pytest.raises(OSError, match="injected publication failure"):
        qp_io.save_result(changed, path, True, common_format=True)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


def test_qp_successful_overwrite_uses_immutable_artifact(tmp_path, qp_result):
    path = tmp_path/"result.json"
    qp_io.save_result(qp_result, path, True)
    old = json.loads(path.read_text())["state_amplitudes"]
    qp_io.save_result(replace(qp_result, energies=qp_result.energies+.1), path, True)
    new = json.loads(path.read_text())["state_amplitudes"]
    assert old["file"] != new["file"]
    assert _sha(tmp_path/old["file"]) == old["sha256"]
    assert _sha(tmp_path/new["file"]) == new["sha256"]
    with np.load(tmp_path/old["file"], allow_pickle=False) as arrays:
        assert_array_equal(arrays["energies"], qp_result.energies)


@pytest.mark.parametrize("backend", ["qp", pytest.param("dmrg", marks=pytest.mark.dmrg)])
def test_unified_persistence_roundtrip_uses_common_records(tmp_path, request, backend):
    result = request.getfixturevalue("qp_result" if backend == "qp" else "mps_result")
    path = tmp_path/"result.json"
    save_result(result, path, save_states=True)
    record = json.loads(path.read_text())
    assert record["format"] == "qdjj-eigenstates" and record["backend"] == backend
    restored = load_result(path)
    assert restored.backend == result.backend == backend
    assert restored.hamiltonian.fingerprint() == result.hamiltonian.fingerprint()
    assert_array_equal(restored.energies, result.energies)
    assert_array_equal(restored.residuals, result.residuals)
    for name, values in result.observables.items():
        assert_allclose(restored.observables[name], values, rtol=0, atol=1e-14)
    if backend == "qp":
        assert_array_equal(restored.vectors, result.vectors)
        assert_array_equal(restored.basis.occupations(), result.basis.occupations())
        assert_array_equal(restored.qp_weights, result.qp_weights)
    else:
        assert_allclose(abs(restored.overlaps(result)), np.eye(2), rtol=0, atol=1e-12)
    save_result(result, path)
    with pytest.raises(ValueError, match="no supported"):
        load_result(path)


def test_qp_loader_preserves_legacy_records_and_four_array_artifacts(tmp_path, qp_result):
    path = tmp_path/"legacy.json"
    qp_io.save_result(qp_result, path, True)
    assert_array_equal(load_result(path).vectors, qp_result.vectors)
    record = json.loads(path.read_text())
    artifact = record["state_amplitudes"]
    state_path = tmp_path/artifact["file"]
    with np.load(state_path, allow_pickle=False) as arrays:
        old_arrays = {name: arrays[name] for name in arrays.files if name != "result_sha256"}
    np.savez_compressed(state_path, **old_arrays)
    artifact["sha256"] = _sha(state_path)
    io.write_json(path, record)
    restored = qp_io.load_result(path)
    assert_array_equal(restored.vectors, qp_result.vectors)
    assert restored.metadata == io.json_value(qp_result.metadata)


@pytest.mark.parametrize("field,value,message", [
    ("file", "../outside.npz", "filename must be local"),
    ("sha256", "0"*64, "checksum mismatch"),
    ("kind", "other", "no supported Fock"),
    ("energies", [100., 101.], "scalar result do not match"),
])
def test_qp_loader_rejects_mismatched_records(tmp_path, qp_result, field, value, message):
    path = tmp_path/"result.json"
    save_result(qp_result, path, True)
    record = json.loads(path.read_text())
    (record if field == "energies" else record["state_artifact"])[field] = value
    io.write_json(path, record)
    with pytest.raises(ValueError, match=message):
        load_result(path)


@pytest.mark.parametrize("fault,message", [
    ("missing", "incomplete Fock"), ("energies", "energies do not match"),
    ("occupations", "do not match the basis"), ("shape", "invalid Fock state vectors"),
    ("nonfinite", "invalid Fock state vectors"), ("normalization", "not orthonormal"),
])
def test_qp_loader_validates_checksummed_numerical_arrays(tmp_path, qp_result, fault, message):
    path = tmp_path/"result.json"
    save_result(qp_result, path, True)
    record = json.loads(path.read_text())
    artifact = record["state_artifact"]
    state_path = tmp_path/artifact["file"]
    with np.load(state_path, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    if fault == "missing":
        del arrays["qp_counts"]
    elif fault == "energies":
        arrays["energies"] += 1.
    elif fault == "occupations":
        arrays["occupations"] = arrays["occupations"][::-1]
    elif fault == "shape":
        arrays["vectors"] = arrays["vectors"][:1]
    elif fault == "nonfinite":
        arrays["vectors"][0, 0] = np.nan
    else:
        arrays["vectors"] *= 2
    np.savez_compressed(state_path, **arrays)
    artifact["sha256"] = record["state_amplitudes"]["sha256"] = _sha(state_path)
    io.write_json(path, record)
    with pytest.raises(ValueError, match=message):
        load_result(path)


@pytest.mark.parametrize("backend", ["qp", pytest.param("dmrg", marks=pytest.mark.dmrg)])
def test_interrupt_after_publication_preserves_referenced_artifact(tmp_path, monkeypatch, request, backend):
    result = request.getfixturevalue("qp_result" if backend == "qp" else "mps_result")
    module = qp_io
    if backend == "dmrg":
        from qdjj_solver.dmrg_solver import io as module
    path = tmp_path/"result.json"
    replace_file = io.os.replace

    def publish_then_interrupt(source, destination):
        replace_file(source, destination)
        raise KeyboardInterrupt

    monkeypatch.setattr(io.os, "replace", publish_then_interrupt)
    with pytest.raises(KeyboardInterrupt):
        module.save_result(result, path, True)
    record = json.loads(path.read_text())
    artifact = record["state_amplitudes" if backend == "qp" else "state_artifact"]
    assert _sha(tmp_path/artifact["file"]) == artifact["sha256"]
    if backend == "qp":
        with np.load(tmp_path/artifact["file"], allow_pickle=False) as arrays:
            assert_array_equal(arrays["energies"], result.energies)
    else:
        assert_array_equal(module.load_result(path).energies, result.energies)


@pytest.mark.dmrg
@pytest.mark.parametrize("layer, value, message", [
    *[("record", value, "unsupported.*result format") for value in (None, [], "result", 1, True)],
    *[("state_artifact", value, "no supported MPS checkpoint") for value in (None, [], "artifact", 1, True)],
    *[("file", value, "checkpoint filename must be a nonempty string") for value in (None, [], {}, 1, True, "")],
])
def test_mps_loaders_and_resume_reject_invalid_record_shapes(tmp_path, capsys, monkeypatch, layer, value, message):
    pytest.importorskip("tenpy")
    from qdjj_solver import cli
    from qdjj_solver.dmrg_solver import io as mps_io

    record = dict(format="qdjj-eigenstates", format_version=1, backend="dmrg",
                  state_artifact=dict(kind="tenpy-mps-hdf5", format_version=1))
    if layer == "record":
        record = value
    elif layer == "file":
        record["state_artifact"]["file"] = value
    else:
        record[layer] = value
    path, config, output = tmp_path/"resume.json", tmp_path/"model.json", tmp_path/"result.json"
    io.write_json(path, record)
    for loader in (load_result, mps_io.load_result):
        with pytest.raises(ValueError, match=message) as error:
            loader(path)
    io.write_json(config, dict(backend="dmrg", model=dict(kind="hamiltonian", **_paired_hamiltonian().record())))
    solver = Mock(side_effect=AssertionError("invalid checkpoint reached the solver"))
    monkeypatch.setattr(cli, "solve", solver)
    assert cli.main(["solve", str(config), "--output", str(output), "--resume", str(path)]) == 2
    captured = capsys.readouterr()
    assert captured.err == f"Calculation failed: {error.value}\n"
    assert not output.exists()
    solver.assert_not_called()


@pytest.mark.dmrg
@pytest.mark.parametrize("suffix", [".pkl", ".pklz", ".mps.h5"])
def test_mps_loader_rejects_pickle_despite_declared_hdf5_format(tmp_path, monkeypatch, suffix):
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    from tenpy.tools import hdf5_io

    # Only inert built-in data; never construct or execute a malicious payload.
    checkpoint = tmp_path/f"checkpoint{suffix}"
    payload = pickle.dumps({"format": "qdjj-mps", "format_version": 1})
    checkpoint.write_bytes(gzip.compress(payload) if suffix == ".pklz" else payload)
    path = tmp_path/"result.json"
    io.write_json(path, dict(format="qdjj-eigenstates", format_version=1, backend="dmrg",
                            state_artifact=dict(kind="tenpy-mps-hdf5", format_version=1,
                                                file=checkpoint.name, sha256=_sha(checkpoint))))
    unpickle = Mock(side_effect=AssertionError("checkpoint selected pickle dispatch"))
    monkeypatch.setattr(pickle, "load", unpickle)
    deserialize = Mock(side_effect=AssertionError("non-HDF5 data reached the TeNPy deserializer"))
    monkeypatch.setattr(hdf5_io, "load_from_hdf5", deserialize)
    with pytest.raises(OSError):
        load_result(path)
    unpickle.assert_not_called()
    deserialize.assert_not_called()


@pytest.mark.dmrg
def test_mps_complex_multistate_roundtrip(tmp_path, mps_result):
    from qdjj_solver.dmrg_solver.io import load_result, save_result

    result = mps_result
    path = tmp_path/"nested"/"result.json"
    save_result(result, path, True)
    record = json.loads(path.read_text())
    artifact = record["state_artifact"]
    assert artifact["kind"] == "tenpy-mps-hdf5" and artifact["format_version"] == 1
    assert artifact["sha256"] == _sha(path.parent/artifact["file"])
    assert Path(artifact["file"]).name == artifact["file"]
    loaded = load_result(path)
    assert_allclose(loaded.energies, [.9-np.sqrt(.94), .9+np.sqrt(.94)], rtol=0, atol=2e-13)
    assert_allclose(abs(loaded.overlaps(result)), np.eye(2), rtol=0, atol=1e-12)
    assert loaded.prepared.coordinate_id == result.prepared.coordinate_id
    assert loaded.hamiltonian.fingerprint() == result.hamiltonian.fingerprint()
    assert loaded.metadata == io.json_value(result.metadata)
    assert loaded.timings == result.timings
    assert_array_equal(io.complex_array(loaded.hamiltonian.metadata["phase"]), [1+2j])
    assert np.max(abs(result.observables["pair"].imag)) > .1
    for name in result.observables:
        assert_allclose(loaded.observables[name], result.observables[name], rtol=0, atol=1e-14)
        operator = result.hamiltonian.observables[name]
        assert_allclose(loaded.matrix_elements(operator), result.matrix_elements(operator), rtol=0, atol=1e-13)
    assert_array_equal(loaded.residuals, result.residuals)


@pytest.mark.dmrg
def test_mps_checkpoint_without_lanczos_probability_tolerance_loads_default(tmp_path, mps_result):
    from qdjj_solver.dmrg_solver import SolverOptions

    metadata = deepcopy(mps_result.metadata)
    del metadata["options"]["lanczos_probability_tolerance"]
    # Save an old-style options record with its matching checkpoint checksum.
    path = tmp_path/"legacy.json"
    save_result(replace(mps_result, metadata=metadata), path, save_states=True)
    restored = load_result(path)
    assert "lanczos_probability_tolerance" not in restored.metadata["options"]
    options = SolverOptions(**restored.metadata["options"])
    assert options.lanczos_probability_tolerance == 1e-14
    assert options.mode_order == list(mps_result.prepared.order)
    assert_allclose(abs(restored.overlaps(mps_result)), np.eye(2), atol=1e-12, rtol=0.)
    assert_array_equal(restored.energies, mps_result.energies)
    assert_array_equal(restored.residuals, mps_result.residuals)


@pytest.mark.dmrg
def test_mps_scalar_only_and_filename_validation(tmp_path, mps_result):
    from qdjj_solver.dmrg_solver.io import load_result, save_result

    with pytest.raises(ValueError, match="end in .json"):
        save_result(mps_result, tmp_path/"bad.h5", True)
    path = tmp_path/"result.json"
    save_result(mps_result, path)
    assert list(tmp_path.iterdir()) == [path]
    with pytest.raises(ValueError, match="no supported MPS"):
        load_result(path)
    path.write_text("{broken")
    with pytest.raises(json.JSONDecodeError):
        load_result(path)
    with pytest.raises(FileNotFoundError):
        load_result(tmp_path/"missing.json")


@pytest.mark.dmrg
@pytest.mark.parametrize("layer, value, message", [
    ("format_version", 2, "unsupported DMRG result"),
    ("kind", "other", "no supported MPS"),
    ("file", "../outside.h5", "filename must be local"),
    ("sha256", "0"*64, "checksum mismatch"),
    ("energies", [100., 101.], "scalar result do not match"),
])
def test_mps_rejects_scalar_and_artifact_corruption(tmp_path, mps_result, layer, value, message):
    from qdjj_solver.dmrg_solver.io import load_result, save_result

    path = tmp_path/"result.json"
    save_result(mps_result, path, True)
    record = json.loads(path.read_text())
    target = record["state_artifact"] if layer in ("kind", "file", "sha256") else record
    target[layer] = value
    io.write_json(path, record)
    with pytest.raises(ValueError, match=message):
        load_result(path)


@pytest.mark.dmrg
@pytest.mark.parametrize("layer, message", [
    ("format", "unsupported MPS checkpoint format"),
    ("hamiltonian", "Hamiltonian fingerprint mismatch"),
    ("coordinate_id", "coordinates do not match"),
    ("states", "state count does not match"),
])
def test_mps_layered_checkpoint_validation(tmp_path, mps_result, layer, message):
    from tenpy.tools import hdf5_io
    from qdjj_solver.dmrg_solver.io import load_result, save_result

    path = tmp_path/"result.json"
    save_result(mps_result, path, True)
    record = json.loads(path.read_text())
    artifact = record.pop("state_artifact")
    checkpoint = tmp_path/artifact["file"]
    data = hdf5_io.load(str(checkpoint))
    if layer == "hamiltonian":
        record["hamiltonian"]["metadata"]["changed"] = True
        data["result_sha256"] = io.fingerprint(record)
    elif layer == "states":
        data["states"] = data["states"][:1]
    else:
        data[layer] = "incorrect"
    hdf5_io.save(data, str(checkpoint))
    artifact["sha256"] = _sha(checkpoint)
    io.write_json(path, dict(record, state_artifact=artifact))
    with pytest.raises(ValueError, match=message):
        load_result(path)


@pytest.mark.dmrg
def test_mps_failed_overwrite_preserves_old_pair(tmp_path, mps_result, monkeypatch):
    from qdjj_solver.dmrg_solver import io as mps_io

    path = tmp_path/"result.json"
    mps_io.save_result(mps_result, path, True)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    changed = replace(mps_result, energies=mps_result.energies+.1)

    def fail(*args, **kwargs):
        raise OSError("injected publication failure")

    # Fail after the replacement HDF5 was written, at JSON publication.
    monkeypatch.setattr(mps_io, "write_json", fail)
    with pytest.raises(OSError, match="injected publication failure"):
        mps_io.save_result(changed, path, True)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
    assert_allclose(mps_io.load_result(path).energies, mps_result.energies, rtol=0, atol=1e-14)


@pytest.mark.dmrg
def test_mps_successful_overwrite_uses_immutable_artifact(tmp_path, mps_result):
    from qdjj_solver.dmrg_solver.io import load_result, save_result

    path = tmp_path/"result.json"
    save_result(mps_result, path, True)
    old_record = json.loads(path.read_text())
    old_artifact = old_record["state_artifact"]
    changed = replace(mps_result, residuals=None)
    save_result(changed, path, True)
    new_artifact = json.loads(path.read_text())["state_artifact"]
    assert old_artifact["file"] != new_artifact["file"]
    assert _sha(tmp_path/old_artifact["file"]) == old_artifact["sha256"]
    assert load_result(path).residuals is None
    # A reader which already saw the old JSON must still be able to restore it.
    snapshot = tmp_path/"snapshot.json"
    io.write_json(snapshot, old_record)
    assert_allclose(load_result(snapshot).residuals, mps_result.residuals, rtol=0, atol=1e-14)
