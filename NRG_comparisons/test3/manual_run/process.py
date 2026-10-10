#!/usr/bin/env python3
"""Standalone, explicit pilot preparation/execution and checked z averaging."""

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
import subprocess

if __package__:
    from .nrg import extract_result, read_param
    from .physics import OBSERVABLES, PRIMARY
else:
    from nrg import extract_result, read_param
    from physics import OBSERVABLES, PRIMARY


HERE = Path(__file__).resolve().parent
ARTIFACTS = ("param", "data", "raw.h5", "DONE", "nrginit.log", "nrg.log")


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def executable(variable):
    value = os.environ.get(variable, "")
    path = Path(value)
    if not path.is_absolute() or not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"set {variable} to an explicit absolute executable path")
    return str(path)  # Preserve symlinks: initializer resources can depend on argv[0].


@contextmanager
def lock(root):
    import fcntl

    with (root / ".manual.lock").open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError(f"another manual workflow owns {root}") from None
        yield


def numbered(root):
    return {path.name for path in root.iterdir() if path.name.isdecimal()}


def common_settings(deck):
    return {"param": {key: value for key, value in deck["param"].items() if key != "z"},
            "extra": deck["extra"]}


def prepare(root):
    """Run only looper with a validated generation-only template, never a solver."""
    root = Path(root).absolute()
    with lock(root):
        if numbered(root) or (root / "sweep.json").exists():
            raise ValueError("existing numbered entries or sweep.json; use a fresh folder")
        text = (root / "param.loop").read_text()
        lines = text.splitlines()
        controls = [line for line in lines if re.match(r"^#(?:AUTOLOOP|LOOP|BEFORE|AFTER|OVERWRITE|CONTINUE|STOP|PRELUDE)", line)]
        prelude = [line for line in controls if re.fullmatch(r"#PRELUDE: \$Nz = [1-9][0-9]*;", line)]
        if len(prelude) != 1 or controls != [prelude[0], "#AUTOLOOP: true"]:
            raise ValueError("param.loop requires literal Nz and preparation-only '#AUTOLOOP: true'")
        nz = int(re.search(r"[0-9]+", prelude[0])[0])
        active = [line for line in lines if line.strip() and not line.startswith("#")]
        if active[:2] != ["@ $i=1; $i<=$Nz; $i++", "!($z=$i/$Nz)"]:
            raise ValueError("param.loop requires the fixed i=1..Nz, z=i/Nz loop")
        # looper evaluates Perl, including its quoted output lines. Allow only
        # literal deck characters after substituting the two approved variables.
        body = "\n".join(active[2:]).replace("$Nz", str(nz)).replace("$z", "1")
        if re.search(r"[^A-Za-z0-9_.,=+\-\[\] \n]", body):
            raise ValueError("param.loop contains unsupported executable expressions")
        subprocess.run([executable("LOOPER"), "param.loop"], cwd=root, check=True)
        if numbered(root) != {str(i) for i in range(1, nz + 1)}:
            raise ValueError("looper did not generate exactly the expected numbered directories")
        entries, common = [], None
        for i in range(1, nz + 1):
            directory = root / str(i)
            if directory.is_symlink() or not directory.is_dir() or (directory / "param").is_symlink():
                raise ValueError("expected ordinary numbered directories and param files")
            if {path.name for path in directory.iterdir()} != {"param"}:
                raise ValueError("looper must generate only param files")
            deck = read_param(directory / "param")
            if int(deck["param"].get("Nz", "0")) != nz or not math.isclose(float(deck["param"]["z"]), i / nz, rel_tol=0, abs_tol=1e-14):
                raise ValueError(f"wrong Nz or z in {directory}/param")
            settings = common_settings(deck)
            if common is not None and common != settings:
                raise ValueError("mixed sweep: all parameters except z must agree")
            common = settings
            entries.append({"directory": str(i), "z": i / nz, "param_sha256": sha256(directory / "param")})
        manifest = {"status": "PILOT ONLY, NOT CONVERGED", "parameters": common,
                    "template_sha256": sha256(root / "param.loop"), "runs": entries}
        write_json(root / "sweep.json", manifest)
        print(f"Prepared {nz} pilot inputs; no initializer or solver launched.", flush=True)
        return manifest


def load_sweep(root):
    root = Path(root)
    manifest = json.loads((root / "sweep.json").read_text())
    nz = int(manifest["parameters"]["param"]["Nz"])
    expected = [str(i) for i in range(1, nz + 1)]
    if ([row["directory"] for row in manifest["runs"]] != expected or numbered(root) != set(expected)
            or sha256(root / "param.loop") != manifest["template_sha256"]):
        raise ValueError("partial/mixed sweep or template changed; use a fresh folder")
    for i, row in enumerate(manifest["runs"], 1):
        directory = root / str(i)
        if directory.is_symlink() or (directory / "param").is_symlink():
            raise ValueError("numbered directories and param files must not be symlinks")
        deck = read_param(directory / "param")
        if (sha256(directory / "param") != row["param_sha256"] or common_settings(deck) != manifest["parameters"]
                or not math.isclose(float(deck["param"]["z"]), i / nz, rel_tol=0, abs_tol=1e-14) or row["z"] != i / nz):
            raise ValueError(f"prepared parameters changed in directory {i}; use a fresh sweep")
    return manifest


def run_one(directory):
    """Explicit execution only; failures remain intact and cannot be resumed."""
    directory = Path(directory).absolute()
    with lock(directory.parent):
        manifest = load_sweep(directory.parent)
        entry = next((row for row in manifest["runs"] if row["directory"] == directory.name), None)
        if entry is None:
            raise ValueError("run_one requires a prepared numbered directory")
        tools = {name: executable(variable) for name, variable in (("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH"))}
        threads = int(os.environ.get("NRG_THREADS", "1"))
        if threads < 1:
            raise ValueError("NRG_THREADS must be a positive integer")
        if {path.name for path in directory.iterdir()} != {"param"}:
            raise ValueError("prior outputs exist; archive this attempt and prepare a fresh sweep")
        started = {"started_utc": datetime.now(timezone.utc).isoformat(), "param_sha256": entry["param_sha256"],
                   "executables": tools, "threads": threads}
        write_json(directory / "started.json", started)
        scratch = directory / "scratch"
        scratch.mkdir()
        environment = dict(os.environ, NRG_WORKDIR=str(scratch), TMPDIR=str(scratch),
                           OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1",
                           MKL_DYNAMIC="FALSE", MKL_THREADING_LAYER="INTEL")
        with (directory / "nrginit.log").open("w") as log:
            subprocess.run([tools["nrginit"], tools["math"]], cwd=directory, env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        if not (directory / "data").is_file() or "Success!" not in (directory / "nrginit.log").read_text():
            raise ValueError("initializer did not produce data and Success!; see nrginit.log")
        with (directory / "nrg.log").open("w") as log:
            subprocess.run([tools["nrg"]], cwd=directory, env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        if any(not (directory / name).is_file() or (directory / name).is_symlink() for name in ARTIFACTS):
            raise ValueError("missing ordinary DONE/raw.h5 or other required outputs")
        if sha256(directory / "param") != entry["param_sha256"]:
            raise ValueError("param changed during execution")
        write_json(directory / "completed.json", dict(started, finished_utc=datetime.now(timezone.utc).isoformat(),
                   artifacts={name: sha256(directory / name) for name in ARTIFACTS}))
        print(f"Completed pilot z={entry['z']:g}; this receipt is not a convergence certificate.", flush=True)


def average(root):
    """Reject partial/mixed sweeps; average only eligible primary state IDs."""
    root = Path(root)
    with lock(root):
        manifest = load_sweep(root)
        records = []
        for entry in manifest["runs"]:
            directory = root / entry["directory"]
            receipt = json.loads((directory / "completed.json").read_text())
            if receipt["param_sha256"] != entry["param_sha256"] or set(receipt["artifacts"]) != set(ARTIFACTS):
                raise ValueError(f"{directory}: incomplete/mismatched completion receipt")
            for name, digest in receipt["artifacts"].items():
                if (directory / name).is_symlink() or sha256(directory / name) != digest:
                    raise ValueError(f"{directory}: output changed since completion: {name}")
            records.append(dict(entry, completion=receipt))
        # Check the whole sweep before spending time on any HDF5 extraction.
        for record in records:
            directory = root / record["directory"]
            result = extract_result(directory)
            primary = set(PRIMARY.values())
            subgap = {row["id"] for row in result["states"] if row["subgap"]}
            result["coverage"]["missing_primary_subgap_ids"] = sorted(primary - subgap)
            result["coverage"]["extra_subgap_ids"] = sorted(subgap - primary)
            record.update(result)
        averages, excluded = [], {}
        for state_id in PRIMARY.values():
            selected, reasons = [], []
            for record in records:
                state = next((row for row in record["states"] if row["id"] == state_id), None)
                reason = "missing" if state is None else "not subgap" if not state["subgap"] else "degenerate" if state["degenerate"] else None
                if reason:
                    reasons.append({"z": record["z"], "reason": reason})
                else:
                    if set(state["values"]) != set(OBSERVABLES) or not all(math.isfinite(v) for v in state["values"].values()):
                        raise ValueError(f"incomplete/nonfinite observables for {state_id}")
                    selected.append(state["values"])
            if reasons:
                excluded[state_id] = reasons
                continue
            for observable in OBSERVABLES:
                values = [row[observable] for row in selected]
                averages.append({"id": state_id, "observable": observable, "nz": len(values),
                                 "mean": math.fsum(values) / len(values), "min": min(values), "max": max(values)})
        summary = {"status": "PILOT ONLY, NOT CONVERGED", "parameters": manifest["parameters"],
                   "averaging": "equal-weight complete z grid; eligible primary IDs only",
                   "continuum_target_established": False, "uncertainty_kind": "not_established",
                   "reference_comparison": None, "per_z": records, "averages": averages,
                   "excluded_primary_ids": excluded,
                   "expected_subgap_coverage_every_twist": all(not row["coverage"]["missing_primary_subgap_ids"]
                                                              and not row["coverage"]["extra_subgap_ids"] for row in records),
                   "limitations": ["No truncation, chain-length, Lambda or Nz convergence established.",
                                   "A single-z pilot is not a z-averaged reference.",
                                   "Extra/missing subgap states are listed per twist, never averaged away.",
                                   "Same-sector degenerate expectations depend on the eigenvector basis.",
                                   "Min/max describe twist variation, not error bars; no transverse observables."]}
        output = root / "processed"
        output.mkdir(exist_ok=True)
        write_json(output / "summary.json", summary)
        with (output / "per_z.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["directory", "z", "id", "twice_sz", "parity", "root",
                                                       "excitation", "subgap", "degenerate", "residual", *OBSERVABLES])
            writer.writeheader()
            for record in records:
                for row in record["states"]:
                    writer.writerow(dict(directory=record["directory"], z=record["z"],
                                         **{k: v for k, v in row.items() if k != "values"}, **row["values"]))
        with (output / "averages.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["id", "observable", "nz", "mean", "min", "max"])
            writer.writeheader()
            writer.writerows(averages)
        print(f"Wrote {output}; excluded primary IDs: {sorted(excluded)}; "
              f"expected subgap coverage: {summary['expected_subgap_coverage_every_twist']}. PILOT ONLY.")
        return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="prepare only by default; execution needs --execute")
    mode = run.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only", action="store_true")
    mode.add_argument("--execute", action="store_true", help="explicitly run the fresh pilot grid sequentially")
    subparsers.add_parser("run-one", help="explicitly run the current prepared numbered directory")
    subparsers.add_parser("average", help="validate every receipt, then extract and average eligible IDs")
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            if args.execute:
                for variable in ("NRGINIT", "NRG", "MATH"):
                    executable(variable)
            manifest = prepare(HERE)
            if args.execute:
                for entry in manifest["runs"]:
                    run_one(HERE / entry["directory"])
        elif args.command == "run-one":
            run_one(Path.cwd())
        else:
            average(HERE)
    except (OSError, ValueError, KeyError, ImportError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"manual test3 NRG: {error}\n")


if __name__ == "__main__":
    main()
