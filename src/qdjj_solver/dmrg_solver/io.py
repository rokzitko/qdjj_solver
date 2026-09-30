"""Portable scalar records and separately checksummed native MPS checkpoints."""

import hashlib
from importlib import import_module
import json
from pathlib import Path

import numpy as np
from tenpy.tools import hdf5_io

from ..common.io import complex_array, fingerprint, result_record, unpublished_artifact, write_json
from ..common.problem import Hamiltonian, Sector
from .mpo import prepare
from .solver import Eigenstates, SolverOptions


def save_result(result, path, save_states=False):
    path = Path(path)
    if path.suffix != ".json":
        raise ValueError("result filename must end in .json")
    record = result_record(result, backend="dmrg")
    if save_states:
        # TeNPy permits scalar-only use without h5py, but its serializer would
        # otherwise raise NameError rather than an actionable import failure.
        import_module("h5py")
        with unpublished_artifact(path, ".mps.h5") as state_path:
            hdf5_io.save(dict(format="qdjj-mps", format_version=1,
                             coordinate_id=result.prepared.coordinate_id,
                             result_sha256=fingerprint(record), states=result.states), str(state_path))
            with state_path.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            record["state_artifact"] = dict(kind="tenpy-mps-hdf5", format_version=1,
                                            file=state_path.name, sha256=digest)
            write_json(path, record)
        return
    write_json(path, record)


def load_result(path):
    """Restore a trusted MPS checkpoint, verifying content and basis identity."""
    path = Path(path)
    record = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(record, dict)
            or (record.get("format"), record.get("format_version"), record.get("backend")) != ("qdjj-eigenstates", 1, "dmrg")):
        raise ValueError("unsupported DMRG result format")
    artifact = record.pop("state_artifact", None)
    if not isinstance(artifact, dict) or artifact.get("kind") != "tenpy-mps-hdf5" or artifact.get("format_version") != 1:
        raise ValueError("this result has no supported MPS checkpoint")
    filename = artifact.get("file")
    if not isinstance(filename, str) or not filename:
        raise ValueError("checkpoint filename must be a nonempty string")
    if Path(filename).name != filename:
        raise ValueError("checkpoint filename must be local to the result directory")
    state_path = path.parent/filename
    with state_path.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != artifact["sha256"]:
            raise ValueError("MPS checkpoint checksum mismatch")
    h5py = import_module("h5py")
    # Enforce the declared format, not filename dispatch. TeNPy data must still be trusted.
    with h5py.File(state_path, "r") as handle:
        data = hdf5_io.load_from_hdf5(handle)
    if data.get("format") != "qdjj-mps" or data.get("format_version") != 1:
        raise ValueError("unsupported MPS checkpoint format")
    if data["result_sha256"] != fingerprint(record):
        raise ValueError("MPS checkpoint and scalar result do not match")
    h = Hamiltonian.from_record(record["hamiltonian"])
    if h.fingerprint() != record["hamiltonian_sha256"]:
        raise ValueError("Hamiltonian fingerprint mismatch")
    calculation = record["calculation"]
    options = SolverOptions(**calculation["options"])
    prepared = prepare(h, Sector(**calculation["sector"]), group_size=options.group_size,
                       mode_order=options.mode_order)
    if prepared.coordinate_id != data["coordinate_id"]:
        raise ValueError("checkpoint coordinates do not match the Hamiltonian and layout")
    energies = np.asarray(record["energies"])
    if len(data["states"]) != len(energies):
        raise ValueError("checkpoint state count does not match the result")
    residuals = None if record["residuals"] is None else np.asarray(record["residuals"])
    return Eigenstates(energies, data["states"], residuals,
                       {name: np.real_if_close(complex_array(values))
                        for name, values in record["observables"].items()},
                       record["timings_seconds"], calculation, h, prepared)
