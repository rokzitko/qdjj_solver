#!/usr/bin/env python3
"""Standalone preparation, checked execution and z averaging for test2."""

import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

from extract import extract_result, read_param


HERE = Path(__file__).resolve().parent
OBSERVABLES = ("signed_gap", "singlet.energy", "doublet.energy", "singlet.P0",
               "singlet.P1", "singlet.P2", "doublet.P0", "doublet.P1",
               "doublet.P2", "doublet.moment")
ARTIFACTS = ("param", "data", "raw.h5", "DONE", "nrginit.log", "nrg.log")


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def executable(name, variable):
    command = os.environ.get(variable, name)
    path = shutil.which(command)
    if path is None:
        raise ValueError(f"missing {name}; load your NRG environment or set {variable} to its executable path")
    # Keep symlinks: some installations use the invoked path to find resources.
    return os.path.abspath(path)


@contextmanager
def lock(path):
    import fcntl

    with path.open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError(f"another manual workflow owns {path}") from None
        yield


def numbered(root):
    return {p.name for p in root.iterdir() if p.name.isdecimal()}


def prepare(root):
    """Use looper only to write inputs, then verify its output independently."""
    with lock(root / ".manual.lock"):
        if numbered(root) or (root / "sweep.json").exists():
            raise ValueError("existing numbered entries or sweep.json; use a fresh folder, not an overwrite/resume sweep")
        template = (root / "param.loop").read_text()
        match = re.search(r"^#PRELUDE:\s*\$Nz\s*=\s*([1-9][0-9]*)\s*;\s*$", template, re.M)
        if match is None:
            raise ValueError("param.loop needs a literal '#PRELUDE: $Nz = 8;' (edit the positive integer)")
        nz = int(match[1])
        directives = [line for line in template.splitlines() if re.match(
            r"#(?:AUTOLOOP|LOOP|BEFORE|AFTER|OVERWRITE|CONTINUE|STOP|PRELUDE)\b", line)]
        if len(directives) != 2 or "#AUTOLOOP: true" not in directives:
            raise ValueError("param.loop must use preparation-only '#AUTOLOOP: true', without other control directives")
        subprocess.run([executable("looper", "LOOPER"), "param.loop"], cwd=root, check=True)
        if numbered(root) != {str(i) for i in range(1, nz + 1)}:
            raise ValueError("looper did not generate exactly the expected numbered directories")
        entries, common = [], None
        for i in range(1, nz + 1):
            directory = root / str(i)
            if directory.is_symlink() or not directory.is_dir() or (directory / "param").is_symlink():
                raise ValueError(f"expected an ordinary directory and param file: {directory}")
            params = read_param(directory / "param")
            if int(params.get("Nz", "0")) != nz or not math.isclose(float(params["z"]), i / nz, rel_tol=0, abs_tol=1e-14):
                raise ValueError(f"wrong Nz or z in {directory}/param")
            settings = {k: v for k, v in params.items() if k != "z"}
            if common is not None and settings != common:
                raise ValueError("all parameters except z must agree across the sweep")
            common = settings
            entries.append({"directory": str(i), "z": i / nz, "param_sha256": sha256(directory / "param")})
        manifest = {"parameters": common, "template_sha256": sha256(root / "param.loop"), "runs": entries}
        write_json(root / "sweep.json", manifest)
        print(f"Prepared {nz} inputs at Lambda={common['Lambda']}; no solver launched.", flush=True)
        return manifest


def load_sweep(root):
    manifest = json.loads((root / "sweep.json").read_text())
    nz = int(manifest["parameters"]["Nz"])
    expected = [str(i) for i in range(1, nz + 1)]
    if ([row["directory"] for row in manifest["runs"]] != expected
            or numbered(root) != set(expected)):
        raise ValueError("incomplete or unexpected numbered directories in the prepared sweep")
    for i, row in enumerate(manifest["runs"], 1):
        directory = root / row["directory"]
        if directory.is_symlink() or (directory / "param").is_symlink():
            raise ValueError("numbered directories and param files must not be symlinks")
        params = read_param(directory / "param")
        if (sha256(directory / "param") != row["param_sha256"]
                or {k: v for k, v in params.items() if k != "z"} != manifest["parameters"]
                or not math.isclose(float(params["z"]), i / nz, rel_tol=0, abs_tol=1e-14)
                or row["z"] != i / nz):
            raise ValueError(f"prepared parameters changed in directory {i}; use a fresh sweep")
    return manifest


def run_one(directory):
    """Run only in a fresh prepared directory; retain failures for inspection."""
    directory = directory.absolute()
    manifest = load_sweep(directory.parent)
    entry = next((row for row in manifest["runs"] if row["directory"] == directory.name), None)
    if entry is None:
        raise ValueError("run_one must be called inside a prepared numbered directory")
    tools = {name: executable(name, variable) for name, variable in
             (("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH"))}
    threads = int(os.environ.get("NRG_THREADS", "4"))
    if threads < 1:
        raise ValueError("NRG_THREADS must be a positive integer")
    with lock(directory / ".run.lock"):
        if any(p.name not in {"param", ".run.lock"} for p in directory.iterdir()):
            raise ValueError(f"{directory}: prior outputs exist; archive this attempt before retrying")
        started = {"started_utc": datetime.now(timezone.utc).isoformat(), "param_sha256": entry["param_sha256"],
                   "executables": tools, "threads": threads}
        write_json(directory / "started.json", started)
        scratch = directory / "scratch"
        scratch.mkdir()
        environment = dict(os.environ, NRG_WORKDIR=str(scratch), TMPDIR=str(scratch),
                           MKL_THREADING_LAYER="INTEL", MKL_NUM_THREADS=str(threads),
                           OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1", MKL_DYNAMIC="FALSE")
        print(f"Running {directory.name}: z={entry['z']:g}", flush=True)
        with (directory / "nrginit.log").open("w") as log:
            subprocess.run([tools["nrginit"], tools["math"]], cwd=directory, env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        if not (directory / "data").is_file() or "Success!" not in (directory / "nrginit.log").read_text():
            raise ValueError(f"{directory}: initializer did not produce data and Success!; see nrginit.log")
        with (directory / "nrg.log").open("w") as log:
            subprocess.run([tools["nrg"]], cwd=directory, env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        if any(not (directory / name).is_file() for name in ARTIFACTS):
            raise ValueError(f"{directory}: missing DONE/raw.h5 or other required output; see nrg.log")
        if sha256(directory / "param") != entry["param_sha256"]:
            raise ValueError("param changed during execution")
        artifacts = {name: sha256(directory / name) for name in ARTIFACTS}
        if (directory / "annotated.dat").exists():
            artifacts["annotated.dat"] = sha256(directory / "annotated.dat")
        write_json(directory / "completed.json", dict(started, artifacts=artifacts,
                   finished_utc=datetime.now(timezone.utc).isoformat()))
        print(f"Completed {directory.name}; raw outputs retained.", flush=True)


def average(root):
    manifest = load_sweep(root)
    reference = json.loads((root / "reference.json").read_text())
    records = []
    for entry in manifest["runs"]:
        directory = root / entry["directory"]
        receipt = json.loads((directory / "completed.json").read_text())
        if (receipt["param_sha256"] != entry["param_sha256"]
                or not set(ARTIFACTS) <= receipt["artifacts"].keys()):
            raise ValueError(f"{directory}: incomplete/mismatched completion receipt")
        for name, digest in receipt["artifacts"].items():
            if name not in (*ARTIFACTS, "annotated.dat") or sha256(directory / name) != digest:
                raise ValueError(f"{directory}: output changed since completion: {name}")
        result = extract_result(directory)
        if set(result["values"]) != set(OBSERVABLES) or not all(math.isfinite(v) for v in result["values"].values()):
            raise ValueError(f"{directory}: incomplete or nonfinite observables")
        records.append(dict(entry, **result, completion=receipt))

    params = manifest["parameters"]
    nz = len(records)
    averages, table = [], []
    # The even members of an eight-twist sweep reproduce the report's Nz=4 grid.
    sizes = [nz] + ([4] if nz > 4 and nz % 4 == 0 else [])
    for count in sizes:
        selected = records[nz // count - 1::nz // count]
        values = {key: sum(row["values"][key] for row in selected) / count for key in OBSERVABLES}
        settings = dict(backend="nrg", nz=count, nmax=int(params["Nmax"]), keep=int(params["keep"]),
                        keepmin=int(params["keepmin"]), keepenergy=float(params["keepenergy"]),
                        **{"lambda": float(params["Lambda"])})
        ref_id = next((key for key, row in reference["nrg_means"].items() if row["settings"] == settings), None)
        averages.append({"settings": settings, "directories": [row["directory"] for row in selected],
                         "values": values, "reference_id": ref_id})
        print(f"\nNz={count}, Lambda={params['Lambda']}; reference: {ref_id or 'no matching stored point'}")
        print(f"{'observable':<20} {'mean':>18} {'report':>18} {'difference':>14}")
        for key, value in values.items():
            ref = reference["nrg_means"][ref_id]["values"][key] if ref_id else None
            difference = value - ref if ref is not None else None
            tolerance = reference["requested_tolerance"]
            ratio = abs(difference) / (tolerance["absolute"] + tolerance["relative"] * abs(value)) if ref is not None else None
            table.append(dict(nz=count, observable=key, mean=value, report_id=ref_id, report=ref,
                              difference=difference, difference_over_requested_tolerance=ratio))
            comparison = f"{ref:18.12g} {difference:14.6g}" if ref is not None else "                 -              -"
            print(f"{key:<20} {value:18.12g} {comparison}")

    summary = dict(parameters=params, averaging="arithmetic mean in increasing z order; complete grids only",
                   uncertainty_kind="not_established", empirical_uncertainties=None,
                   continuum_target_established=False, reference_source=reference["source"],
                   per_z=records, averages=averages, comparisons=table)
    output = root / "processed"
    output.mkdir(exist_ok=True)
    write_json(output / "summary.json", summary)
    with (output / "per_z.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["directory", "z", *OBSERVABLES])
        writer.writeheader()
        writer.writerows(dict(directory=row["directory"], z=row["z"], **row["values"]) for row in records)
    with (output / "averages.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    print(f"\nWrote {output}; differences/requested tolerance are NOT error estimates.")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="prepare and run a fresh sweep, sequentially")
    run.add_argument("--prepare-only", action="store_true", help="write numbered param files without launching NRG")
    subparsers.add_parser("run-one", help="run in the current prepared numbered directory")
    subparsers.add_parser("average", help="validate every output and print/write averages")
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            if not args.prepare_only:
                for name, variable in (("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH")):
                    executable(name, variable)
            manifest = prepare(HERE)
            if not args.prepare_only:
                for entry in manifest["runs"]:
                    run_one(HERE / entry["directory"])
        elif args.command == "run-one":
            run_one(Path.cwd())
        else:
            average(HERE)
    except (OSError, ValueError, KeyError, ImportError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"manual NRG: {error}\n")


if __name__ == "__main__":
    main()
