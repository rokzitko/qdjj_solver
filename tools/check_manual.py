"""Check a built manual: TeX diagnostics, readable PDF, identity and version.

    cd docs/manual && make pdf && cd ../.. && python tools/check_manual.py

Reads ``docs/manual/build/manual.log`` and ``manual.blg``. Underfull boxes and the
font-shape and microtype notices of the class are ignored; everything else that LaTeX or
biber reports as a warning fails the check, as do undefined references and citations.
"""

import argparse
from pathlib import Path
import re
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT/"docs"/"manual"/"build"
IGNORED = ("Font shape", "microtype", "Size substitutions")


def problems(build=BUILD, basename="manual"):
    log = (build/f"{basename}.log").read_text(errors="replace")
    found = []
    for line in log.splitlines():
        if re.search(r"(LaTeX|Package \w+|Class \w+) Warning", line) and not any(i in line for i in IGNORED):
            found.append(line.strip())
        elif line.startswith("Overfull"):
            found.append(line.strip())
        elif line.startswith("Missing character:"):
            found.append(line.strip())
        elif "undefined" in line and ("Reference" in line or "Citation" in line):
            found.append(line.strip())
    pages = re.search(r"Output written on \S+ \((\d+) pages?", log)
    if not pages:
        found.append("no PDF written")
    blg = build/f"{basename}.blg"
    if blg.is_file():
        found += [line.strip() for line in blg.read_text(errors="replace").splitlines()
                  if re.search(r"\b(WARN|ERROR)\b", line)]
    return found, int(pages.group(1)) if pages else 0


def artifact_problems(path, *, title=None, version=None, pages=None, required=()):
    """Inspect the actual artifact; a successful-looking log is not sufficient."""
    from pypdf import PdfReader
    if not path.is_file() or not path.stat().st_size:
        return [f"missing or empty PDF: {path.name}"]
    try:
        reader = PdfReader(path, strict=True)
        count = len(reader.pages)
        if count == 0:
            return ["PDF has no pages"]
        found = []
        if pages is not None and count != pages:
            found.append(f"PDF has {count} pages but log reports {pages}")
        identity = str((reader.metadata or {}).get("/Title", ""))
        if title and title.casefold() not in identity.casefold():
            found.append(f"unexpected PDF title: {identity}")
        front = "\n".join(page.extract_text() or "" for page in reader.pages[:4])
        if version and version not in front:
            found.append(f"software version {version} absent from front matter")
        if required:
            text = " ".join(" ".join((page.extract_text() or "").split()) for page in reader.pages)
            for phrase in required:
                if phrase not in text:
                    found.append(f"required PDF content absent: {phrase}")
        return found
    except Exception as error:
        return [f"unreadable PDF: {type(error).__name__}: {error}"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--build-dir", type=Path, default=BUILD)
    parser.add_argument("--basename", default="manual")
    parser.add_argument("--title", default="physics manual")
    parser.add_argument("--version", default=tomllib.loads((ROOT/"pyproject.toml").read_text())["project"]["version"])
    args = parser.parse_args(argv)
    found, pages = problems(args.build_dir, args.basename)
    found += artifact_problems(args.build_dir/f"{args.basename}.pdf", title=args.title, version=args.version, pages=pages)
    for problem in found:
        print(problem)
    print(f"{pages} pages; {len(found)} problem(s)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
