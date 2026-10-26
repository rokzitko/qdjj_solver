"""Run every Python listing of the physics manual, in document order, in one namespace.

    python tools/run_manual_listings.py          # all listings (several minutes)
    python tools/run_manual_listings.py --list   # show the listings without running them

Listings build on earlier ones, as a reader would type them. Lines ``%! code``
placed directly before ``\\begin{lstlisting}`` hold code that the listing needs
but the text does not show (for example a frequency grid); ``%! skip`` excludes
a listing. Listings with ``language=bash`` are not Python and are skipped.
The listings run in a temporary directory with one numerical thread.
"""

import argparse
import os
from pathlib import Path
import re
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
MANUAL = ROOT/"docs"/"manual"
BEGIN = re.compile(r"\\begin\{lstlisting\}(\[[^\]]*\])?")


def chapters(manual=MANUAL):
    """Chapter files in the order of ``\\input`` in manual.tex."""
    text = (manual/"manual.tex").read_text()
    return [manual/f"{name}.tex" for name in re.findall(r"^\\input\{(chapters/[^}]+)\}", text, re.M)]


def listings(manual=MANUAL):
    """[(path, line, code)] of the Python listings, with the code of their ``%!`` lines prepended."""
    found = []
    active = set()
    def parse(path):
        if path in active:
            raise ValueError(f"recursive manual input: {path}")
        active.add(path)
        lines = path.read_text().split("\n")
        index = 0
        while index < len(lines):
            included = re.match(r"\s*\\input\{(generated/[^}]+)\}", lines[index])
            if included:
                source = (manual/(included.group(1)+".tex")).resolve()
                if not source.is_relative_to(manual.resolve()):
                    raise ValueError(f"manual input outside source directory: {source}")
                parse(source)
            match = BEGIN.search(lines[index])
            if not match:
                index += 1
                continue
            options = match.group(1) or ""
            hidden, back = [], index - 1
            while back >= 0 and lines[back].lstrip().startswith("%!"):
                hidden.insert(0, lines[back].lstrip()[2:].strip())
                back -= 1
            end = next(i for i in range(index + 1, len(lines)) if "\\end{lstlisting}" in lines[i])
            body = lines[index + 1:end]
            if "language=bash" not in options and "skip" not in hidden:
                found.append((path, index + 1, "\n".join(hidden + body)))
            index = end + 1
        active.remove(path)
    for path in chapters(manual):
        parse(path)
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--list", action="store_true", help="list the listings without running them")
    arguments = parser.parse_args(argv)
    found = listings()
    if arguments.list:
        for path, line, code in found:
            print(f"{path.relative_to(ROOT)}:{line}: {len(code.splitlines())} lines")
        return 0
    from check_tutorials import THREADS  # the thread variables, shared with the tutorial checks
    for name in THREADS:
        os.environ[name] = "1"
    namespace = {"__name__": "__manual__"}
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        for path, line, code in found:
            where = f"{path.relative_to(ROOT)}:{line}"
            begin = time.perf_counter()
            try:
                exec(compile(code, where, "exec"), namespace)
            except Exception as error:  # report the listing that failed, then stop
                print(f"FAILED {where}: {type(error).__name__}: {error}")
                return 1
            print(f"ok {where} ({time.perf_counter() - begin:.1f} s)")
    print(f"{len(found)} listings ran in {time.perf_counter() - start:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
