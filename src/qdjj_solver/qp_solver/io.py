"""Versioned, human-readable records of Hamiltonians and numerical results."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ..common.io import (bath_from_config as bath_from_config, json_value as json_value,
                         complex_array, configured, fingerprint, model_from_config, numerical_environment as numerical_environment,
                         unpublished_artifact, validate_config, write_json as write_json)
from ..common.problem import Hamiltonian
from .basis import FockBasis, Sector
from .solver import Eigenstates, SolverOptions


def result_record(result, include_hamiltonian=False):
    record = dict(format="bcs-qp-eigenstates", format_version=1,
                energies=result.energies, residuals=result.residuals,
                qp_weights=result.qp_weights, observables=result.observables,
                timings_seconds=result.timings, calculation=result.metadata,
                environment=numerical_environment())
    if include_hamiltonian:
        if result.hamiltonian is None:
            raise ValueError("the result has no Hamiltonian specification")
        record["hamiltonian"] = result.hamiltonian.record()
    return record


def save_result(result, path, save_vectors=False, *, common_format=False):
    """Save JSON; optionally save amplitudes and occupations in a companion NPZ."""
    path = Path(path)
    if path.suffix != ".json":
        raise ValueError("result filename must end in .json")
    if common_format:
        from ..common.io import result_record as common_result_record
        record = common_result_record(result, backend="qp")
        record["qp_weights"] = result.qp_weights
    else:
        record = result_record(result, include_hamiltonian=True)
    if save_vectors:
        with unpublished_artifact(path, ".npz") as vector_path:
            np.savez_compressed(vector_path, energies=result.energies, vectors=result.vectors,
                                occupations=result.basis.occupations(), qp_counts=result.basis.qp_counts(),
                                result_sha256=fingerprint(record))
            with vector_path.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            record["state_amplitudes"] = dict(file=vector_path.name, sha256=digest)
            if common_format:
                record["state_artifact"] = dict(kind="fock-amplitudes-npz", format_version=1,
                                                file=vector_path.name, sha256=digest)
            write_json(path, record)
        return
    write_json(path, record)


def load_result(path):
    """Restore common or legacy QP records, checking the artifact and Fock basis."""
    path = Path(path)
    record = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(record, dict) or record.get("format_version") != 1:
        raise ValueError("unsupported QP result format")
    if record.get("format") == "qdjj-eigenstates" and record.get("backend") == "qp":
        artifact = record.get("state_artifact")
        if (not isinstance(artifact, dict) or artifact.get("kind") != "fock-amplitudes-npz"
                or artifact.get("format_version") != 1):
            raise ValueError("this result has no supported Fock state artifact")
    elif record.get("format") == "bcs-qp-eigenstates":
        artifact = record.get("state_amplitudes")
        if not isinstance(artifact, dict):
            raise ValueError("this result has no supported Fock state artifact")
    else:
        raise ValueError("unsupported QP result format")
    if Path(artifact["file"]).name != artifact["file"]:
        raise ValueError("state filename must be local to the result directory")
    state_path = path.parent/artifact["file"]
    with state_path.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != artifact["sha256"]:
            raise ValueError("Fock state artifact checksum mismatch")
    with np.load(state_path, allow_pickle=False) as arrays:
        required = {"energies", "vectors", "occupations", "qp_counts"}
        if not required <= set(arrays.files):
            raise ValueError("incomplete Fock state artifact")
        # Earlier NPZ files contain the four numerical arrays only. New saves
        # additionally bind the full scalar record, as the MPS serializer does.
        if "result_sha256" in arrays:
            scalar = {key: value for key, value in record.items()
                      if key not in ("state_artifact", "state_amplitudes")}
            if arrays["result_sha256"].shape != () or arrays["result_sha256"].item() != fingerprint(scalar):
                raise ValueError("Fock state artifact and scalar result do not match")
        saved_energies = arrays["energies"]
        vectors, occupations, counts = arrays["vectors"], arrays["occupations"], arrays["qp_counts"]
    h = Hamiltonian.from_record(record["hamiltonian"])
    if record["format"] == "qdjj-eigenstates" and h.fingerprint() != record["hamiltonian_sha256"]:
        raise ValueError("Hamiltonian fingerprint mismatch")
    calculation = record["calculation"]
    sector = Sector(**calculation["sector"])
    h.check_sector(sector)
    options = SolverOptions(**calculation["options"])
    basis = FockBasis(h.nimp, h.spins, calculation["qp_cutoff"], sector, h.eta_labels, options.max_dimension)
    energies = np.asarray(record["energies"])
    if energies.ndim != 1 or len(energies) == 0 or not np.array_equal(energies, saved_energies):
        raise ValueError("Fock state energies do not match the scalar record")
    if (occupations.dtype != np.dtype("uint64") or not np.array_equal(occupations, basis.occupations())
            or not np.array_equal(counts, basis.qp_counts())):
        raise ValueError("Fock state occupations or QP counts do not match the basis")
    if (vectors.shape != (basis.dimension, len(energies)) or vectors.dtype.kind not in "fc"
            or not np.all(np.isfinite(vectors))):
        raise ValueError("invalid Fock state vectors")
    if np.linalg.norm(vectors.conj().T @ vectors-np.eye(len(energies))) > 1e-8:
        raise ValueError("Fock states are not orthonormal")
    residuals = np.asarray(record["residuals"])
    weights = np.asarray(record["qp_weights"])
    observables = {name: np.real_if_close(complex_array(values)) for name, values in record["observables"].items()}
    shapes = [(energies, (len(energies),)), (residuals, energies.shape),
              (weights, (len(energies), basis.cutoff+1))]
    shapes.extend((values, energies.shape) for values in observables.values())
    if any(values.shape != shape or values.dtype.kind not in "fciub" or not np.all(np.isfinite(values))
           for values, shape in shapes):
        raise ValueError("invalid QP scalar result arrays")
    return Eigenstates(energies, vectors, residuals, basis, weights, observables,
                       record["timings_seconds"], calculation, h)


def calculation_from_config(config):
    """Interpret a reference-junction or a general multi-orbital JSON model."""
    config = validate_config(config, allow_backend=False)
    sector = configured(Sector, config.get("sector", {}), "sector")
    options = configured(SolverOptions, config.get("solver", {}), "solver")
    hamiltonian = model_from_config(config)
    return hamiltonian, config.get("cutoff", 2), sector, options
