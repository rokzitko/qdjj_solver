"""Backend-independent configuration, JSON encoding, and installed-code provenance."""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import contextmanager
import hashlib
from importlib import metadata
from inspect import signature
import json
import os
from pathlib import Path
import platform
import tempfile

import numpy as np

from .. import __version__


def json_value(value):
    """Convert NumPy values and complex scalars/arrays to portable JSON values."""
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(v) for v in value]
    return value


def complex_array(value):
    """Read a real array or the explicit ``{real: ..., imag: ...}`` encoding."""
    if isinstance(value, dict):
        if set(value) != {"real", "imag"}:
            raise ValueError("a complex array needs real and imag fields")
        real, imag = np.asarray(value["real"]), np.asarray(value["imag"])
        if real.shape != imag.shape:
            raise ValueError("real and imaginary array shapes differ")
        return real + 1j*imag
    return np.asarray(value, dtype=complex)


def fingerprint(value):
    text = json.dumps(json_value(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_json(path, data):
    """Write atomically, rejecting non-finite values before replacing a record."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(json_value(data), indent=2, sort_keys=True, allow_nan=False) + "\n"
    name = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as handle:
            name = handle.name
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if name is not None and os.path.exists(name):
            os.unlink(name)


@contextmanager
def unpublished_artifact(path, suffix):
    """Reserve a unique companion file for an atomically published JSON record.

    The caller writes the artifact and publishes its reference with ``write_json``
    inside this context. Ordinary failures remove the new file. Interruptions
    retain it because the JSON reference may already have been published.
    Older artifacts remain valid for readers that already obtained the old JSON.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f"{path.stem}.", suffix=suffix, dir=path.parent)
    artifact = Path(name)
    os.close(descriptor)
    try:
        yield artifact
    except Exception:
        artifact.unlink(missing_ok=True)
        raise


def numerical_environment():
    packages = {}
    for name in ("qdjj-solver", "numpy", "scipy", "physics-tenpy", "h5py", "threadpoolctl"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    return dict(qdjj_solver=__version__, qpsolver=__version__,
                python=platform.python_version(), packages=packages,
                numpy=np.__version__, scipy=packages.get("scipy"),
                system=platform.system(), release=platform.release(), machine=platform.machine())


def source_provenance():
    """Hash the implementation actually imported, including an installed native core.

    Relative module paths are portable between editable and wheel installations;
    no checkout location, user home directory, or Git working tree is recorded.
    """
    root = Path(__file__).resolve().parents[1]
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix in (".py", ".so", ".pyd", ".dylib"):
            with path.open("rb") as handle:
                files[path.relative_to(root).as_posix()] = hashlib.file_digest(handle, "sha256").hexdigest()
    compatibility = root.parent/"qpsolver"
    if compatibility.is_dir():
        for path in sorted(compatibility.glob("*.py")):
            with path.open("rb") as handle:
                files[f"compat/qpsolver/{path.name}"] = hashlib.file_digest(handle, "sha256").hexdigest()
    # With an editable install the extension may be in site-packages rather than src/.
    import sys
    core = sys.modules.get("qdjj_solver.qp_solver._core")
    if core is not None:
        path = Path(core.__file__)
        with path.open("rb") as handle:
            files[f"qp_solver/{path.name}"] = hashlib.file_digest(handle, "sha256").hexdigest()
    return dict(version=__version__, sha256=fingerprint(files), files=files)


def _config_object(value, path, *, allowed=None, required=()):
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise ValueError(f"{path} must be an object with string keys")
    result = dict(value)
    if allowed is not None and result.keys()-set(allowed):
        raise ValueError(f"unknown {path} keys: {sorted(result.keys()-set(allowed))}")
    if set(required)-result.keys():
        raise ValueError(f"missing {path} keys: {sorted(set(required)-result.keys())}")
    return result


def _config_array(value, path):
    if (not isinstance(value, (list, tuple, np.ndarray))
            or isinstance(value, np.ndarray) and value.ndim == 0):
        raise ValueError(f"{path} must be an array")
    return list(value)


def configured(constructor, values, path, *args):
    """Validate constructor keywords and retain the configuration path in errors."""
    allowed = list(signature(constructor).parameters)[len(args):]
    values = _config_object(values, path, allowed=allowed)
    try:
        return constructor(*args, **values)
    except (ValueError, TypeError, KeyError) as error:
        raise ValueError(f"{path}: {error}") from error


def validate_config(config, *, allow_backend=True):
    """Validate the versioned envelope without importing either numerical backend."""
    allowed = {"format_version", "model", "bath", "cutoff", "sector", "solver"}
    if allow_backend:
        allowed.add("backend")
    config = _config_object(config, "configuration", allowed=allowed)
    version = config.get("format_version", 1)
    if type(version) is not int or version != 1:
        raise ValueError("unsupported configuration version")
    for name in ("model", "bath", "sector", "solver"):
        if name in config:
            config[name] = _config_object(config[name], name)
    cutoff = config.get("cutoff")
    if cutoff is not None and (not isinstance(cutoff, (int, np.integer))
                               or isinstance(cutoff, bool) or cutoff < 0):
        raise ValueError("cutoff must be a nonnegative integer or null")
    sector = config.get("sector", {})
    for name in ("parity", "twice_sz", "eta", "particle_number"):
        if name in sector and sector[name] is not None:
            if not isinstance(sector[name], (int, np.integer)) or isinstance(sector[name], bool):
                raise ValueError(f"sector.{name} must be an integer")
    for section, names in (("model", ("compress", "symmetry", "eta_basis", "compress_pairs")),
                           ("solver", ("mixer", "calculate_residuals", "require_convergence"))):
        for name in names:
            if name in config.get(section, {}) and not isinstance(config[section][name], bool):
                raise ValueError(f"{section}.{name} must be a boolean")
    return config


def _operator_from_config(records, path):
    from .algebra import FermionOperator
    operator = FermionOperator()
    for index, record in enumerate(_config_array(records, path)):
        location = f"{path}[{index}]"
        record = _config_object(record, location, allowed=("operators", "real", "imag"),
                                required=("operators", "real", "imag"))
        word = _config_array(record["operators"], f"{location}.operators")
        try:
            operator += FermionOperator({tuple(word): complex(record["real"], record["imag"])})
        except (ValueError, TypeError) as error:
            raise ValueError(f"{location}: {error}") from error
    return operator


def _observables_from_config(values, path):
    return {name: _operator_from_config(records, f"{path}.{name}")
            for name, records in _config_object(values, path).items()}


def _contacts_from_config(values, path):
    contacts = {}
    for index, record in enumerate(_config_array(values, path)):
        location = f"{path}[{index}]"
        record = _config_object(record, location, allowed=("leads", "matrix"), required=("leads", "matrix"))
        leads = _config_array(record["leads"], f"{location}.leads")
        if (len(leads) != 2 or any(not isinstance(i, (int, np.integer)) or isinstance(i, bool) for i in leads)
                or not 0 <= leads[0] < leads[1]):
            raise ValueError(f"{location}.leads must contain two integer indices with 0 <= l < m")
        pair = tuple(leads)
        if pair in contacts:
            raise ValueError(f"{location}.leads duplicates contact {pair}")
        try:
            contacts[pair] = complex_array(record["matrix"])
        except (ValueError, TypeError) as error:
            raise ValueError(f"{location}.matrix: {error}") from error
    return contacts


def bath_from_config(config, *, path="bath"):
    from .baths import DiscreteBath, chain_expansion, cosh_grid, fit_surrogate
    config = _config_object(config, path)
    kind = config.pop("kind", "cosh-grid")
    if "metadata" in config:
        config["metadata"] = _config_object(config["metadata"], f"{path}.metadata")
    if kind == "cosh-grid":
        return configured(cosh_grid, config, path)
    if kind == "surrogate":
        return configured(fit_surrogate, config, path)
    if kind == "chain-expansion":
        return configured(chain_expansion, config, path)
    if kind == "discrete":
        return configured(DiscreteBath, config, path)
    raise ValueError(f"{path}: unknown bath representation {kind!r}")


def model_from_config(config, *, default_compress=True):
    """Construct one shared physical model; no numerical backend is imported."""
    from .models import Impurity, make_model, reference_model
    from .problem import Hamiltonian
    config = validate_config(config)
    model = dict(config.get("model", {}))
    kind = model.pop("kind", "reference")
    if kind == "reference":
        bath = bath_from_config(config.get("bath", {"pairs": 8}))
        model.setdefault("compress", default_compress)
        return configured(reference_model, model, "model", bath)
    if kind == "hamiltonian":
        model = _config_object(model, "model", allowed=signature(Hamiltonian).parameters,
                               required=("operator", "nimp", "spins"))
        model["operator"] = _operator_from_config(model["operator"], "model.operator")
        model["spins"] = _config_array(model["spins"], "model.spins")
        if model.get("eta_labels") is not None:
            model["eta_labels"] = _config_array(model["eta_labels"], "model.eta_labels")
        model["observables"] = _observables_from_config(model.get("observables", {}), "model.observables")
        model["metadata"] = _config_object(model.get("metadata", {}), "model.metadata")
        return configured(Hamiltonian, model, "model")
    if kind == "multi-orbital":
        allowed = set(signature(make_model).parameters)-{"baths"} | {"reservoirs"}
        model = _config_object(model, "model", allowed=allowed, required=("impurity", "reservoirs", "tunneling"))
        impurity = _config_object(model.pop("impurity"), "model.impurity",
                                  allowed=("orbitals", "terms", "observables", "metadata"), required=("orbitals", "terms"))
        imp = configured(Impurity, dict(orbitals=impurity["orbitals"],
                         operator=_operator_from_config(impurity["terms"], "model.impurity.terms"),
                         observables=_observables_from_config(impurity.get("observables", {}), "model.impurity.observables"),
                         metadata=_config_object(impurity.get("metadata", {}), "model.impurity.metadata")), "model.impurity")
        baths = [bath_from_config(b, path=f"model.reservoirs[{i}]")
                 for i, b in enumerate(_config_array(model.pop("reservoirs"), "model.reservoirs"))]
        hopping = []
        for i, t in enumerate(_config_array(model.pop("tunneling"), "model.tunneling")):
            try:
                hopping.append(complex_array(t))
            except (ValueError, TypeError) as error:
                raise ValueError(f"model.tunneling[{i}]: {error}") from error
        model["direct"] = _contacts_from_config(model.get("direct", []), "model.direct")
        if "direct_derivatives" in model:
            model["direct_derivatives"] = {
                name: _contacts_from_config(value, f"model.direct_derivatives.{name}")
                for name, value in _config_object(model["direct_derivatives"], "model.direct_derivatives").items()}
        return configured(make_model, model, "model", imp, baths, hopping)
    raise ValueError("model kind must be reference, multi-orbital, or hamiltonian")


def result_record(result, *, backend, include_hamiltonian=True):
    """Common scalar-result envelope; backend-specific state files are separate."""
    record = dict(format="qdjj-eigenstates", format_version=1, backend=backend,
                  energies=result.energies, residuals=result.residuals,
                  observables=result.observables, timings_seconds=result.timings,
                  calculation=result.metadata, environment=numerical_environment(),
                  implementation=source_provenance(),
                  hamiltonian_sha256=result.hamiltonian.fingerprint())
    if include_hamiltonian:
        record["hamiltonian"] = result.hamiltonian.record()
    return record
