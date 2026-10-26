"""Run the tutorials in examples/ and compare their results with the checked references.

    python tools/check_tutorials.py                      # quick profile of every tutorial
    python tools/check_tutorials.py --profile full       # full profiles (minutes; some need extras)
    python tools/check_tutorials.py 08_two_tone          # selected tutorials

Each tutorial ``examples/NN_topic/run.py`` accepts ``--profile quick|full`` and
``--output DIR`` and writes ``DIR/results.json`` with ``{"values": {key: number}}``.
``reference/<profile>.json`` lists the expected values: ``{"value", "atol", "rtol"}``
for a number, or ``{"max": bound}`` for an upper bound such as a residual.
"""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import tomllib

ROOT = Path(__file__).resolve().parents[1]
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def catalog(root=ROOT):
    return tomllib.loads((root/"examples"/"catalog.toml").read_text())["tutorial"]


def compare(reference, results):
    """Problems of ``results`` (the values dict of results.json) against a reference record."""
    problems = []
    for key, expected in reference["values"].items():
        if key not in results:
            problems.append(f"{key}: missing from results")
            continue
        value = results[key]
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            problems.append(f"{key}: {value!r} is not a finite number")
        elif "max" in expected:
            if value > expected["max"]:
                problems.append(f"{key}: {value:.6g} exceeds {expected['max']:.6g}")
        else:
            tolerance = expected.get("atol", 0.) + expected.get("rtol", 0.)*abs(expected["value"])
            if abs(value - expected["value"]) > tolerance:
                problems.append(f"{key}: {value!r} differs from {expected['value']!r} by more than {tolerance:.3g}")
    return problems


def run_tutorial(directory, profile, output, root=ROOT, timeout=None):
    """Run one tutorial profile in a subprocess with one numerical thread; return results.json."""
    environment = dict(os.environ, **{name: "1" for name in THREADS}, PYTHONDONTWRITEBYTECODE="1")
    command = [sys.executable, "-B", str(root/"examples"/directory/"run.py"), "--profile", profile,
               "--output", str(output)]
    completed = subprocess.run(command, cwd=root, env=environment, capture_output=True, text=True, timeout=timeout)
    if completed.returncode:
        raise RuntimeError(f"{directory} {profile} failed:\n{completed.stdout[-2000:]}\n{completed.stderr[-4000:]}")
    return json.loads((Path(output)/"results.json").read_text())


def check(directory, profile, output, root=ROOT, timeout=None):
    reference = json.loads((root/"examples"/directory/"reference"/f"{profile}.json").read_text())
    results = run_tutorial(directory, profile, output, root, timeout)
    if results.get("profile") != profile:
        return [f"results.json reports profile {results.get('profile')!r}"]
    return compare(reference, results["values"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("tutorials", nargs="*", help="tutorial directories (default: all in examples/catalog.toml)")
    parser.add_argument("--profile", choices=("quick", "full"), default="quick", help="profile to run")
    parser.add_argument("--keep", type=Path, help="keep the outputs in this directory")
    arguments = parser.parse_args(argv)
    known = [entry["dir"] for entry in catalog()]
    selected = arguments.tutorials or known
    unknown = sorted(set(selected) - set(known))
    if unknown:
        parser.error(f"not in examples/catalog.toml: {', '.join(unknown)}")
    failed = 0
    with tempfile.TemporaryDirectory() as scratch:
        for directory in selected:
            output = (arguments.keep or Path(scratch))/directory/arguments.profile
            start = time.perf_counter()
            try:
                problems = check(directory, arguments.profile, output)
            except RuntimeError as error:
                problems = [str(error)]
            elapsed = time.perf_counter() - start
            status = "ok" if not problems else "FAILED"
            print(f"{directory} [{arguments.profile}] {status} ({elapsed:.1f} s)")
            for problem in problems:
                print(f"  {problem}")
            failed += bool(problems)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
