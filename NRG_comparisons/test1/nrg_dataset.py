"""Content-verified, raw-replayed NRG inputs, independent of campaign analysis.

Snapshot paths are relative to ``run.HERE / 'runs'``, not to the snapshot file.
Only freeze() audits original archives; load_snapshot() needs the frozen input
and its caller-pinned SHA, even when the original archives are no longer readable.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time

from NRG_comparisons.test1 import nrg, run


NUMERICAL_FIELDS = ("values", "finite_converged", "diagnostics", "finite_convergence_meaning")
RECORD_FIELDS = ("task", "physical", "status", *NUMERICAL_FIELDS)
REQUIRED_ARTIFACTS = {"param", "data", "nrginit.log", "nrg.log", "raw.h5"}


def fingerprint():
    """Return the pure postprocessor's code/runtime identity; resolve no solvers."""
    from qdjj_solver.common.io import numerical_environment

    environment = numerical_environment()
    return dict(
        postprocess_nrg_sha256=hashlib.sha256(inspect.getsource(run.postprocess_nrg).encode()).hexdigest(),
        flatten_sha256=hashlib.sha256(inspect.getsource(run.flatten).encode()).hexdigest(),
        nrg_sha256=hashlib.sha256(Path(nrg.__file__).read_bytes()).hexdigest(),
        environment={k: environment[k] for k in ("python", "numpy", "system", "release", "machine")} |
                    {"packages": {k: environment["packages"][k] for k in ("numpy", "h5py")}})


def _decode(contents):
    value = json.loads(contents)
    # json.loads otherwise accepts NaN and Infinity, which are not audit evidence.
    json.dumps(value, allow_nan=False)
    if not isinstance(value, dict):
        raise ValueError("NRG archive must contain a JSON object")
    return value


def load_snapshot(path, expected_sha256):
    """Read just a pinned snapshot, without provenance discovery or raw replay."""
    contents = Path(path).read_bytes()
    if hashlib.sha256(contents).hexdigest() != expected_sha256:
        raise ValueError("NRG snapshot SHA256 mismatch")
    result = _decode(contents)
    if result.get("format_version") != 1 or len(result["records"]) != len(result["sources"]):
        raise ValueError("NRG snapshot format or record/source alignment mismatch")
    return result


def _delta(saved, actual, where="replay"):
    """Compare all diagnostics, with absolute-only float tolerance and exact scalars."""
    if isinstance(saved, dict) and isinstance(actual, dict) and saved.keys() == actual.keys():
        return max((_delta(saved[k], actual[k], f"{where}.{k}") for k in saved), default=0.)
    if isinstance(saved, list) and isinstance(actual, list) and len(saved) == len(actual):
        return max((_delta(a, b, f"{where}[{i}]") for i, (a, b) in enumerate(zip(saved, actual, strict=True))), default=0.)
    if type(saved) is float and type(actual) is float:
        difference = abs(saved - actual)
        if math.isfinite(saved) and math.isfinite(actual) and difference <= 1e-12:
            return difference
    elif type(saved) is type(actual) and not isinstance(saved, (dict, list)) and saved == actual:
        return 0.
    raise ValueError(f"NRG {where} mismatch")


def freeze(source_campaign, physical, output_path):
    """Audit a transitive archive and atomically publish a *new* JSON input.

    Source campaigns and every referenced file must stay within this case/runs.
    The output parent must already exist; existing outputs and writes inside any
    source campaign are refused. Legacy nrg_reuse and targeted baseline edges are
    followed, with all producer manifests pinned. Each distinct completed NRG
    leaf is content-hashed and replayed, never executed. Different completed
    attempts for the same exact task are ambiguous even if their values agree.

    Records retain all numerical diagnostics but no solver-time accounting.
    sources[i] pins records[i]'s original full result, request, task, producer and
    artifacts (annotated.dat=None records absence). verification_bytes counts
    bytes read for content hashing, excluding additional reads during replay and
    runtime discovery; verification_seconds is audit time, not solver time.
    """
    started = time.monotonic()
    runs = (run.HERE / "runs").resolve()
    output = Path(output_path).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"NRG snapshot output already exists: {output}")
    if not output.parent.is_dir():
        raise FileNotFoundError(f"NRG snapshot output parent must exist: {output.parent}")
    output = output.resolve()

    def inside(path):
        path = Path(path).resolve()
        if not path.is_relative_to(runs):
            raise ValueError("NRG archive paths must stay within this case/runs")
        return path

    def relative(base, path):
        if not isinstance(path, str) or Path(path).is_absolute():
            raise ValueError("NRG archive references must be relative within this case/runs")
        return inside(base / path)

    source = inside(source_campaign)
    hashes, documents = {}, {}
    manifests, snapshots, leaves, visited, active = {}, {}, {}, {}, set()
    verification_bytes = 0
    unpinned = object()

    def content(path, expected=unpinned, *, document=False):
        nonlocal verification_bytes
        if expected is not unpinned and (not isinstance(expected, str) or re.fullmatch(r"[0-9a-f]{64}", expected) is None):
            raise ValueError("NRG reference requires a SHA256 hash")
        path = inside(path)
        if path not in hashes:
            checksum = hashlib.sha256()
            chunks = [] if document else None
            with path.open("rb") as handle:
                while block := handle.read(1024 * 1024):
                    checksum.update(block)
                    verification_bytes += len(block)
                    if chunks is not None:
                        chunks.append(block)
            hashes[path] = checksum.hexdigest()
            if document:
                documents[path] = _decode(b"".join(chunks))
        if expected is not unpinned and hashes[path] != expected:
            raise ValueError(f"NRG SHA256 mismatch: {path}")
        return documents[path] if document else hashes[path]

    def pin(path):
        path = inside(path)
        return dict(path=str(path.relative_to(runs)), sha256=content(path))

    current = run.nrg_runtime(run.provenance(["nrg"], nrg_only=True))
    processor = fingerprint()
    if processor["environment"] != current["environment"] or processor["nrg_sha256"] != current["nrg_sha256"]:
        raise ValueError("NRG postprocessor/runtime identity mismatch")
    identity = dict(format_version=1, energy_unit="Delta", physical=physical)

    def manifest(path, expected=unpinned):
        path = inside(path)
        if path.name != "manifest.json" or path.parent.parent != runs:
            raise ValueError("NRG producer manifest must identify a campaign directly under case/runs")
        if output.is_relative_to(path.parent):
            raise ValueError("NRG snapshot output must not modify a source campaign")
        value = content(path, expected, document=True)
        if any(value["manifest"].get(k) != v for k, v in identity.items()):
            raise ValueError("NRG manifest physical identity or units/version mismatch")
        if run.nrg_runtime(value["provenance"]) != current:
            raise ValueError("NRG producer runtime/renderer provenance mismatch")
        manifests[path] = pin(path)
        return value

    def leaf(path, producer, expected=unpinned):
        path, producer = inside(path), inside(producer)
        record = content(path, expected, document=True)
        owner = manifest(producer)
        task = record["task"]
        key = run.digest(task)
        directory = path.parent
        if (path.name != "result.json" or directory.parent.parent != producer.parent / "tasks"
                or directory.parent.name != "nrg-" + key[:16]
                or re.fullmatch(r"attempt-[0-9]+", directory.name) is None):
            raise ValueError("NRG result path/producer/task identity mismatch")
        if (task.get("backend") != "nrg" or task not in owner["tasks"]
                or record.get("status") != "completed" or record.get("physical") != physical):
            raise ValueError("NRG result task/physical/completion identity mismatch")
        if path in leaves:
            if leaves[path][1]["producer_manifest"] != pin(producer):
                raise ValueError("NRG result has conflicting producers")
            return path
        task_path, request_path = directory.parent / "task.json", directory / "request.json"
        if content(task_path, document=True) != task:
            raise ValueError("NRG task.json mismatch")
        request = content(request_path, document=True)
        if request.get("task") != task or request.get("physical") != physical:
            raise ValueError("NRG request identity mismatch")
        deck = nrg.render_param(physical, task["numerical"], task["units"], task["clean"])
        content(directory / "param", hashlib.sha256(deck.encode()).hexdigest())
        artifacts = record["artifacts"]
        if not REQUIRED_ARTIFACTS <= artifacts.keys():
            raise ValueError("NRG missing required artifact hashes")
        checked = {}
        for name, checksum in artifacts.items():
            if Path(name).name != name or name in (".", ".."):
                raise ValueError("NRG artifact names must be plain filenames")
            checked[name] = content(directory / name, checksum)
        annotated = inside(directory / "annotated.dat")
        checked["annotated.dat"] = content(annotated) if annotated.exists() else None
        saved = {k: record[k] for k in NUMERICAL_FIELDS}
        replay = run.postprocess_nrg(task, physical, directory)
        maximum = _delta(saved, replay)
        normalized = {k: record[k] for k in RECORD_FIELDS}
        evidence = dict(**pin(path), task_sha256=key, producer_manifest=pin(producer),
                        task_file=pin(task_path), request=pin(request_path), artifacts=checked,
                        replay_max_delta=maximum)
        leaves[path] = normalized, evidence
        return path

    def unique(paths):
        by_task = {}
        for path in paths:
            record, _ = leaves[path]
            key = run.digest(record["task"])
            if key in by_task and by_task[key] != path:
                other = leaves[by_task[key]][0]
                kind = "conflicting numerical" if any(record[k] != other[k] for k in NUMERICAL_FIELDS) else "ambiguous"
                raise ValueError(f"NRG {kind} duplicate task attempts: {by_task[key]} and {path}")
            by_task[key] = path
        return list(by_task.values())

    def baseline(path, expected):
        path = inside(path)
        if path in active:
            raise ValueError("NRG archive reference cycle")
        value = content(path, expected, document=True)
        if (value.get("format_version") != 1 or value.get("runs_root") != str(runs)
                or value.get("identity") != identity or len(value["records"]) != len(value["sources"])):
            raise ValueError("NRG baseline snapshot format/root/physical identity mismatch")
        if value["postprocessor"] != processor or value["runtime"] != current:
            raise ValueError("NRG baseline postprocessor/runtime mismatch")
        active.add(path)
        snapshots[path] = pin(path)
        pins = {}
        for reference in value["manifests"]:
            producer = relative(runs, reference["path"])
            manifest(producer, reference["sha256"])
            pins[producer] = reference["sha256"]
        root = value["source_manifest"]
        root_path = relative(runs, root["path"])
        if pins.get(root_path) != root["sha256"]:
            raise ValueError("NRG baseline source manifest is not pinned")
        available = visit(root_path, root["sha256"])
        imported = []
        for record, reference in zip(value["records"], value["sources"], strict=True):
            original = relative(runs, reference["path"])
            producer_ref = reference["producer_manifest"]
            producer = relative(runs, producer_ref["path"])
            if pins.get(producer) != producer_ref["sha256"] or original not in available:
                raise ValueError("NRG baseline leaf/producer is not in the pinned source graph")
            leaf(original, producer, reference["sha256"])
            normalized, evidence = leaves[original]
            if record != normalized or any(reference[k] != evidence[k] for k in evidence if k != "replay_max_delta"):
                raise ValueError("NRG baseline record/evidence mismatch")
            imported.append(original)
        if set(imported) != set(available):
            raise ValueError("NRG baseline records do not cover the pinned source graph")
        # Every transitively visited manifest must have been pinned by this input.
        for producer in graph_manifests[root_path]:
            if pins.get(producer) != hashes[producer]:
                raise ValueError("NRG baseline transitive manifest is not pinned")
        active.remove(path)
        return unique(imported), set(pins)

    graph_manifests = {}

    def visit(path, expected=unpinned):
        path = inside(path)
        if path in active:
            raise ValueError("NRG archive reference cycle")
        value = manifest(path, expected)
        if path in visited:
            return visited[path]
        active.add(path)
        paths, graph = [], {path}
        campaign = path.parent
        if "baseline" in value and "nrg_reuse" in value:
            raise ValueError("NRG campaign cannot mix baseline and legacy reuse edges")
        if "baseline" in value:
            reference = value["baseline"]
            imported, ancestors = baseline(relative(campaign, reference["path"]), reference["sha256"])
            paths.extend(imported)
            graph.update(ancestors)
        if "nrg_reuse" in value:
            reuse = value["nrg_reuse"]
            parent = relative(campaign, reuse["manifest_path"])
            if parent in active:
                raise ValueError("NRG archive reference cycle")
            archived = manifest(parent, reuse["manifest_sha256"])
            if reuse["provenance"] != archived["provenance"]:
                raise ValueError("NRG recorded source provenance differs from actual manifest")
            available = visit(parent, reuse["manifest_sha256"])
            graph.update(graph_manifests[parent])
            for reference in reuse["results"]:
                original = relative(campaign, reference["path"])
                if original not in available:
                    raise ValueError("NRG reused result is not in the referenced producer graph")
                content(original, reference["sha256"])
                record, evidence = leaves[original]
                if reference["task_sha256"] != evidence["task_sha256"] or record["task"] not in value["tasks"]:
                    raise ValueError("NRG reused task identity mismatch")
                paths.append(original)
        for candidate in sorted((campaign / "tasks").glob("*/attempt-*/result.json")):
            result = content(candidate, document=True)
            if result.get("task", {}).get("backend") == "nrg" and result.get("status") == "completed":
                paths.append(leaf(candidate, path))
        active.remove(path)
        visited[path] = unique(paths)
        graph_manifests[path] = graph
        return visited[path]

    paths = visit(source / "manifest.json")
    result = dict(format_version=1, runs_root=str(runs), identity=identity,
                  source_manifest=pin(source / "manifest.json"),
                  records=[leaves[p][0] for p in paths], sources=[leaves[p][1] for p in paths],
                  manifests=[manifests[p] for p in sorted(manifests)],
                  snapshots=[snapshots[p] for p in sorted(snapshots)],
                  postprocessor=processor, runtime=current,
                  verification_bytes=verification_bytes, verification_files=len(hashes),
                  verification_seconds=time.monotonic() - started)
    serialized = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    # Hard-link publication is atomic and refuses a racing existing destination.
    # The temporary file is private, in the output parent, never in an archive.
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, prefix=".nrg-snapshot-") as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())
        os.link(handle.name, output)
    return result
