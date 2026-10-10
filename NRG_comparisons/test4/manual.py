#!/usr/bin/env python3
"""Export a portable test4 NRG bundle; preparation never runs NRG or Mathematica."""

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
import tempfile

if __package__:
    from . import nrg, physics
else:
    import nrg
    import physics


HERE = Path(__file__).resolve().parent
MODELS = ("model.m", "modeloperators.m")
SOURCES = (*MODELS, "nrg.py", "physics.py", "manual.py")
ARTIFACTS = ("param", *MODELS, "data", "raw.h5", "DONE", "nrginit.log", "nrg.log")
NUMERICAL = {"lambda": 4., "z": 1., "nmax": 20, "keep": 200, "keepmin": 40, "keepenergy": 14.}
STATUS = "UNCONVERGED pilot template; no continuum accuracy established"
README = r"""# Manual Test4 NRG

**UNCONVERGED pilot template.** Defaults: Nz=8, Lambda=4, Nmax=20,
keep=200, keepmin=40, keepenergy=14. These are not production settings.
The physical defaults come from the exported `physics.PHYSICAL`.
This bundle contains snapshots of the canonical model, not an independent model.
Copy the whole folder anywhere before preparing; no source-checkout paths are used.

```sh
./1_run                    # prepare only, equivalent to manual.py run
(cd 1 && ../run_one)        # execute one twist; repeat for every numbered directory
./2_proc                   # average only a complete, checked grid
```

`run_one` uses the current prepared numbered directory, not the bundle root.
Alternatively, in a fresh unprepared bundle, `./1_run --execute` prepares and
executes the entire grid sequentially. It does not automatically average.
All wrappers accept `PYTHON=/path/to/python3`; Python 3.11+ is required.
Preparation needs `looper`; execution needs `nrginit`, `math`, and `nrg` on PATH,
or explicit executable overrides `LOOPER`, `NRGINIT`, `MATH`, and `NRG`.
Extraction uses NumPy and h5py, imported only when averaging.
`python -B manual.py run --prepare-only` is also supported.

Edit only `param.loop` before preparation. Keep the preparation-only looper
directives and integer i=1..Nz loop, with z=i/Nz. Physical/numerical deck values
must be literal real numbers, never Mathematica `!` expressions. Lambda >= 4
is mandatory. The paper L/R model supports real dimensional WN and WS together;
other model representations are not supported by this manual workflow.
Root source hashes are checked against `bundle.json`; each numbered directory
gets its own model.m and modeloperators.m for the initializer's CWD guards.
The template, all parameters except the prescribed z, and source hashes are
frozen in `sweep.json`. Do not mix directories or edit models after export.

Execution refuses stale outputs and retries, requires initializer data and
`Success!`, then NRG `DONE` and `raw.h5`, and records artifact hashes in
`completed.json`. Failure leaves the attempt intact; use a fresh bundle.
A root lock serializes workflows in this folder only, not other folders.
Scratch belongs to each numbered directory. `NRG_THREADS` defaults to 1 and
must be 1..4 (OpenMP/MKL settings; OpenBLAS is set to 1).
There is **no wall-time or memory limit**, no global Mathematica license lock,
and no cleanup or killing of unrelated processes. Budget resources explicitly.

## Units And Outputs

Default deck energies are in positive gap units, Delta=1; bandrescale=100 is
the normal half-bandwidth D, Gamma=0.4 is total hybridization (0.2 per lead).
The single `bcsgap=-Delta/D` sets **both** physical L/R pairing coefficients
negative. Do not replace it with opposite signs or separate bcsgap1/2 inputs.
nu and phi are dimensionless (phi in radians); U, Bz, WN and WS are energies,
not rho W. keepenergy is in Wilson-shell units.

`processed/per_z.csv` retains every extracted row, including all final-shell
subgap roots and the first three roots per saved sector. Extra subgap roots
and same-sector degeneracies are flagged, not discarded. `processed/means.csv`
contains equal-weight means only for primary IDs (even, odd_down, odd_up)
present, subgap, and nondegenerate at **every** grid point. `summary.json`
retains all per-z states, coverage, excluded IDs, receipts and provenance.
All 14 observables are included. Missing runs or changed hashes abort averaging;
an existing `processed/` is not overwritten.

`energy` is centered and `energy_paper` uses the paper impurity constant, both
in units of Delta after subtracting the two disconnected centered lead vacua:

```math
E_{\mathrm{paper}}/\Delta=E_{\mathrm{centered}}/\Delta+U(\nu-1)^2/(2\Delta).
```

`current` is the real expectation of dH/dphi divided by Delta, namely
$I/(2e\Delta/\hbar)$, with positive sign for phi=phi_R-phi_L. Values in
$e\Delta/\hbar$ units would be twice as large. Spin observables are longitudinal,
in hbar units. Twist min/max are not error bars. Retention, chain length,
Lambda, Nz and continuum convergence remain unestablished; complete saved
subgap coverage is not proof of continuum completeness.
"""


def sha256(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"missing ordinary file: {path}")
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    with path.open("x") as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_hashes(root):
    expected = json.loads((root / "bundle.json").read_text())["source_sha256"]
    actual = {name: sha256(root / name) for name in SOURCES}
    if actual != expected:
        raise ValueError("exported source/model hashes changed; export a fresh bundle")
    return actual


def export(destination=None):
    """Copy canonical sources only at export time; never overwrite a target."""
    destination = Path(destination) if destination is not None else HERE / "manual_run"
    destination = destination.absolute()
    if destination.is_symlink() or destination.exists() or HERE.is_relative_to(destination.resolve()):
        raise ValueError("export requires a fresh destination, not canonical source or an existing target")
    sources = {name: sha256(HERE / name) for name in SOURCES}
    if (HERE / "bundle.json").exists():
        source_hashes(HERE)
    deck = nrg.render_param(physics.PHYSICAL, NUMERICAL)
    deck = deck.replace("[param]\n", "[param]\nNz=$Nz\n")
    deck = re.sub(r"^z=.*$", "z=$z", deck, flags=re.M)
    template = (f"#!/usr/bin/env looper\n# {STATUS}\n#PRELUDE: $Nz = 8;\n"
                "#AUTOLOOP: true\n@ $i=1; $i<=$Nz; $i++\n!($z=$i/$Nz)\n" + deck)
    destination.mkdir()
    for name in SOURCES:
        shutil.copy2(HERE / name, destination / name)
        if sha256(destination / name) != sources[name]:
            raise ValueError(f"source changed during export: {name}")
    (destination / "param.loop").write_text(template)
    (destination / "README.md").write_text(README)
    for name, command in (("1_run", "run"), ("run_one", "run-one"), ("2_proc", "average")):
        wrapper = destination / name
        wrapper.write_text('#!/bin/sh\nset -eu\nexec "${PYTHON:-python3}" -B '
                           f'"$(dirname "$0")/manual.py" {command} "$@"\n')
        wrapper.chmod(0o755)
    write_json(destination / "bundle.json", {"status": STATUS, "source_sha256": sources})
    print(f"Exported {destination}; no calculations performed.", flush=True)
    return destination


def executable(name, variable):
    path = shutil.which(os.environ.get(variable, name))
    if path is None:
        raise ValueError(f"missing {name}; set {variable} or load it on PATH")
    # Preserve symlinks: initializer resources can depend on the invoked path.
    return os.path.abspath(path)


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


def prepare(root, *, execute=False):
    """Prepare only by default; hold the root lock across an explicit full run."""
    root = Path(root).resolve()
    with lock(root):
        sources = source_hashes(root)
        if numbered(root) or any((root / name).exists() for name in ("sweep.json", "processed")):
            raise ValueError("existing sweep/numbered entries; use a fresh folder")
        template_hash = sha256(root / "param.loop")
        lines = (root / "param.loop").read_text().splitlines()
        controls = [line for line in lines if re.match(
            r"^#(?:AUTOLOOP|LOOP|BEFORE|AFTER|OVERWRITE|CONTINUE|STOP|PRELUDE)", line)]
        if (len(controls) != 2 or not re.fullmatch(r"#PRELUDE: \$Nz = [1-9][0-9]*;", controls[0])
                or controls[1] != "#AUTOLOOP: true"):
            raise ValueError("require literal Nz PRELUDE and preparation-only '#AUTOLOOP: true'")
        nz = int(re.search(r"[0-9]+", controls[0])[0])
        active = [line for line in lines if line.strip() and not line.startswith("#")]
        if (active[:2] != ["@ $i=1; $i<=$Nz; $i++", "!($z=$i/$Nz)"]
                or active.count("Nz=$Nz") != 1 or active.count("z=$z") != 1):
            raise ValueError("require integer i=1..Nz, z=i/Nz loop and Nz=$Nz, z=$z deck entries")
        body = "\n".join(active[2:]).replace("Nz=$Nz", f"Nz={nz}").replace("z=$z", "z=1") + "\n"
        # Looper interpolates Perl strings. Reject executable syntax BEFORE launch,
        # then let the canonical adapter validate all literal physics/numerics.
        if re.search(r"[^A-Za-z0-9_.,=+\-\[\] \n]", body):
            raise ValueError("param.loop requires literal values, not executable expressions")
        with tempfile.TemporaryDirectory(prefix=".validate-", dir=root) as temporary:
            path = Path(temporary) / "param"
            path.write_text(body)
            expected = nrg.read_param(path)
        if expected["representation"] != "paper":
            raise ValueError("manual test4 supports only the canonical paper model")
        tools = execution_tools() if execute else None
        subprocess.run([executable("looper", "LOOPER"), "param.loop"], cwd=root, check=True)
        if numbered(root) != {str(i) for i in range(1, nz + 1)}:
            raise ValueError("looper did not generate exactly the complete numbered grid")
        entries = []
        for i in range(1, nz + 1):
            directory = root / str(i)
            if directory.is_symlink() or not directory.is_dir() or {p.name for p in directory.iterdir()} != {"param"}:
                raise ValueError("looper must generate ordinary directories containing only param")
            param_hash = sha256(directory / "param")
            deck = nrg.read_param(directory / "param")
            if (common_settings(deck) != common_settings(expected)
                    or not math.isclose(deck["numerical"]["z"], i / nz, rel_tol=0, abs_tol=1e-14)):
                raise ValueError("mixed settings or incorrect z grid from looper")
            for name in MODELS:
                shutil.copy2(root / name, directory / name)
                if sha256(directory / name) != sources[name]:
                    raise ValueError(f"model changed during preparation: {name}")
            entries.append({"directory": str(i), "z": i / nz, "param_sha256": param_hash,
                            "physical": deck["physical"], "numerical": deck["numerical"]})
        if source_hashes(root) != sources or sha256(root / "param.loop") != template_hash:
            raise ValueError("sources/template changed during preparation")
        manifest = {"status": STATUS, "parameters": common_settings(expected), "source_sha256": sources,
                    "template_sha256": template_hash, "runs": entries}
        write_json(root / "sweep.json", manifest)
        print(f"Prepared {nz} inputs; no initializer or solver launched during preparation.", flush=True)
        if execute:
            for entry in entries:
                _run_one(root / entry["directory"], tools)
        return manifest


def load_sweep(root):
    sources = source_hashes(root)
    manifest = json.loads((root / "sweep.json").read_text())
    nz = int(manifest["parameters"]["param"]["Nz"])
    expected = [str(i) for i in range(1, nz + 1)]
    if (nz < 1 or [row["directory"] for row in manifest["runs"]] != expected
            or numbered(root) != set(expected) or manifest["source_sha256"] != sources
            or sha256(root / "param.loop") != manifest["template_sha256"]):
        raise ValueError("partial/mixed sweep or changed template/sources")
    for i, entry in enumerate(manifest["runs"], 1):
        directory = root / entry["directory"]
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError("expected an ordinary numbered directory")
        if sha256(directory / "param") != entry["param_sha256"]:
            raise ValueError(f"prepared param changed: {directory}")
        deck = nrg.read_param(directory / "param")
        if (common_settings(deck) != manifest["parameters"] or entry["physical"] != deck["physical"]
                or entry["numerical"] != deck["numerical"] or entry["z"] != i / nz
                or not math.isclose(deck["numerical"]["z"], i / nz, rel_tol=0, abs_tol=1e-14)):
            raise ValueError(f"mixed physical/numerical inputs: {directory}")
        for name in MODELS:
            if sha256(directory / name) != sources[name]:
                raise ValueError(f"prepared model changed: {directory / name}")
    return manifest


def execution_tools():
    threads = int(os.environ.get("NRG_THREADS", "1"))
    if not 1 <= threads <= 4:
        raise ValueError("NRG_THREADS must be an integer in 1..4")
    return {"executables": {name: executable(name, variable) for name, variable in
                            (("nrginit", "NRGINIT"), ("math", "MATH"), ("nrg", "NRG"))}, "threads": threads}


def _run_one(directory, tools):
    """Caller holds the bundle lock, including between sequential grid points."""
    manifest = load_sweep(directory.parent)
    entry = next((row for row in manifest["runs"] if row["directory"] == directory.name), None)
    if entry is None:
        raise ValueError("run-one requires the current prepared numbered directory")
    if {path.name for path in directory.iterdir()} != {"param", *MODELS}:
        raise ValueError("prior outputs exist; archive this attempt and use a fresh sweep")
    started = dict(tools, started_utc=datetime.now(timezone.utc).isoformat(),
                   param_sha256=entry["param_sha256"], source_sha256=manifest["source_sha256"],
                   template_sha256=manifest["template_sha256"])
    write_json(directory / "started.json", started)
    scratch = directory / "scratch"
    scratch.mkdir()
    environment = dict(os.environ, NRG_WORKDIR=str(scratch), TMPDIR=str(scratch),
                       OMP_NUM_THREADS=str(tools["threads"]), MKL_NUM_THREADS=str(tools["threads"]),
                       OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_DYNAMIC="FALSE",
                       OMP_DYNAMIC="FALSE", MKL_THREADING_LAYER="INTEL")
    commands = tools["executables"]
    with (directory / "nrginit.log").open("x") as log:
        subprocess.run([commands["nrginit"], commands["math"]], cwd=directory, env=environment,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    sha256(directory / "data")
    if "Success!" not in (directory / "nrginit.log").read_text():
        raise ValueError("initializer missing Success!; see nrginit.log")
    if load_sweep(directory.parent) != manifest:
        raise ValueError("inputs changed during initialization")
    with (directory / "nrg.log").open("x") as log:
        subprocess.run([commands["nrg"]], cwd=directory, env=environment,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    artifacts = {name: sha256(directory / name) for name in ARTIFACTS}
    if load_sweep(directory.parent) != manifest:
        raise ValueError("inputs changed during execution")
    write_json(directory / "completed.json", dict(started, artifacts=artifacts,
               finished_utc=datetime.now(timezone.utc).isoformat()))
    print(f"Completed z={entry['z']:g}; receipt is not a convergence certificate.", flush=True)


def run_one(directory):
    directory = Path(directory).absolute()
    with lock(directory.parent):
        _run_one(directory, execution_tools())


def average(root):
    """Check every receipt before extraction; never average a primary subgrid."""
    root = Path(root).resolve()
    with lock(root):
        manifest = load_sweep(root)
        records = []
        for entry in manifest["runs"]:
            directory = root / entry["directory"]
            receipt = json.loads((directory / "completed.json").read_text())
            if (receipt["param_sha256"] != entry["param_sha256"]
                    or receipt["source_sha256"] != manifest["source_sha256"]
                    or receipt["template_sha256"] != manifest["template_sha256"]
                    or set(receipt["artifacts"]) != set(ARTIFACTS)):
                raise ValueError(f"incomplete/mismatched completion receipt: {directory}")
            for name, digest in receipt["artifacts"].items():
                if sha256(directory / name) != digest:
                    raise ValueError(f"output changed since completion: {directory / name}")
            if "Success!" not in (directory / "nrginit.log").read_text():
                raise ValueError(f"initializer missing Success!: {directory}")
            records.append(dict(entry, completion=receipt))
        primary = set(physics.PRIMARY.values())
        for record in records:
            result = nrg.extract_result(root / record["directory"])
            subgap = {row["id"] for row in result["states"] if row["subgap"]}
            result["coverage"].update(missing_primary_subgap_ids=sorted(primary - subgap),
                                      extra_subgap_ids=sorted(subgap - primary))
            # Keep the frozen deck inputs too: extraction reports physical gap units.
            record["result"] = result
        means, excluded = [], {}
        for state_id in physics.PRIMARY.values():
            selected, reasons = [], []
            for record in records:
                state = next((row for row in record["result"]["states"] if row["id"] == state_id), None)
                reason = ("missing" if state is None else "not subgap" if not state["subgap"]
                          else "degenerate" if state["degenerate"] else None)
                if reason:
                    reasons.append({"z": record["z"], "reason": reason})
                else:
                    selected.append(state["values"])
            if reasons:
                excluded[state_id] = reasons
                continue
            for observable in physics.OBSERVABLES:
                values = [row[observable] for row in selected]
                means.append({"id": state_id, "observable": observable, "nz": len(values),
                              "mean": math.fsum(values) / len(values), "min": min(values), "max": max(values)})
        summary = {"status": STATUS, "sweep": manifest, "per_z": records, "means": means,
                   "excluded_primary_ids": excluded, "continuum_completeness_established": False,
                   "averaging": "equal-weight complete z grid; only primary, subgap, nondegenerate IDs",
                   "twist_min_max_are_error_bars": False,
                   "expected_three_subgap_states_every_twist": all(
                       row["result"]["coverage"]["expected_three_subgap_states"] for row in records)}
        output = root / "processed"
        output.mkdir()
        write_json(output / "summary.json", summary)
        with (output / "per_z.csv").open("x", newline="") as stream:
            fields = ["directory", "z", "id", "twice_sz", "parity", "root", "excitation", "subgap",
                      "degenerate", "residual", "primary", "extra_subgap", *physics.OBSERVABLES]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for record in records:
                for row in record["result"]["states"]:
                    writer.writerow(dict(directory=record["directory"], z=record["z"],
                                         **{k: v for k, v in row.items() if k != "values"}, **row["values"],
                                         primary=row["id"] in primary,
                                         extra_subgap=row["subgap"] and row["id"] not in primary))
        with (output / "means.csv").open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["id", "observable", "nz", "mean", "min", "max"])
            writer.writeheader()
            writer.writerows(means)
        print(f"Wrote {output}; excluded primary IDs: {sorted(excluded)}. {STATUS}", flush=True)
        return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    export_parser = commands.add_parser("export", help="export a fresh standalone bundle; no calculations")
    export_parser.add_argument("destination", nargs="?", type=Path, default=HERE / "manual_run")
    run = commands.add_parser("run", help="prepare only by default")
    mode = run.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only", action="store_true")
    mode.add_argument("--execute", action="store_true", help="explicitly execute a fresh complete grid sequentially")
    commands.add_parser("run-one", help="execute the current prepared numbered directory")
    commands.add_parser("average", help="check all receipts and extract the complete grid")
    args = parser.parse_args(argv)
    try:
        if args.command == "export":
            export(args.destination)
        elif args.command == "run":
            prepare(HERE, execute=args.execute)
        elif args.command == "run-one":
            run_one(Path.cwd())
        else:
            average(HERE)
    except (OSError, ValueError, KeyError, ImportError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"manual test4 NRG: {error}\n")


if __name__ == "__main__":
    main()
