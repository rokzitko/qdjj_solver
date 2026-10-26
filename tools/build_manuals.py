"""Build checked, versioned PDF manuals from the current source tree and saved data.

python -B tools/build_manuals.py --manual all
python -B tools/build_manuals.py --manual user

The docs, docs-pdf, plots and dmrg extras and the documented system tools are required.
Release builds additionally specify the exact tag, source commit and source archive.
"""

import argparse
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tomllib

from check_manual import artifact_problems, problems

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def source_identity(root=ROOT):
    if not (root/".git").exists():
        return None, False
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    dirty = bool(subprocess.check_output(["git", "-C", str(root), "status", "--porcelain"], text=True).strip())
    return commit, dirty


def verify_versions(tag=None):
    import qdjj_solver
    source = tomllib.loads((ROOT/"pyproject.toml").read_text())["project"]["version"]
    installed = metadata.version("qdjj-solver")
    if not source == installed == qdjj_solver.__version__:
        raise ValueError(f"source/runtime/installed version mismatch: {source}, {qdjj_solver.__version__}, {installed}")
    if tag is not None and tag != "v"+source:
        raise ValueError(f"release tag {tag} does not match version {source}")
    return source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manual", choices=("user", "physics", "all"), default="all")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"build"/"manuals")
    parser.add_argument("--release-tag")
    parser.add_argument("--source-ref", help="full source commit for release provenance")
    parser.add_argument("--source-archive", type=Path)
    args = parser.parse_args(argv)
    if args.release_tag and (args.manual != "all" or not args.source_ref or not args.source_archive):
        parser.error("release builds require --manual all, --source-ref and --source-archive")
    if args.source_archive and not args.source_archive.is_file():
        parser.error("--source-archive must name the archive used for this build")
    version = verify_versions(args.release_tag)
    commit, dirty = source_identity()
    commit = args.source_ref or commit
    if args.release_tag and dirty:
        parser.error("release documentation requires an unmodified source tree")
    label = (f"Release {args.release_tag}" if args.release_tag else
             f"Development snapshot {(commit or 'unversioned source')[:12]}"+(" (local changes)" if dirty else ""))
    environment = dict(os.environ, QDJJ_DOCS_REF=args.release_tag or commit or f"v{version}", QDJJ_DOCS_LABEL=label)
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[name] = "1"
    output = args.output_dir.resolve()
    logs = ROOT/"build"/"manual-logs"
    logs.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)

    def run(name, command):
        print(f"{name}: {' '.join(command)}", flush=True)
        with (logs/f"{name}.log").open("w") as handle:
            completed = subprocess.run(command, cwd=ROOT, env=environment, stdout=handle,
                                       stderr=subprocess.STDOUT, timeout=600)
        if completed.returncode:
            raise RuntimeError(f"{name} failed; see {logs/f'{name}.log'}")

    run("generated", [sys.executable, "-B", "tools/gen_docs.py", "--check"])
    run("release-version", [sys.executable, "-B", "tools/check_release.py"]+
        (["--tag", args.release_tag] if args.release_tag else []))
    selected = ("physics", "user") if args.manual == "all" else (args.manual,)
    artifacts = {}
    for manual in selected:
        if manual == "physics":
            run("physics", ["make", "-C", "docs/manual", "pdf", f"PYTHON={shlex.quote(sys.executable)} -B"])
            folder, basename, identity, suffix = ROOT/"docs/manual/build", "manual", "physics manual", "physics-manual"
            required = ("From a Hamiltonian to a checked calculation",)
        else:
            run("user-sphinx", [sys.executable, "-m", "sphinx", "-W", "--keep-going", "-t", "handbook",
                                "-c", "docs", "-b", "latex", ".", "build/user-manual"])
            run("user-tex", ["make", "-C", "build/user-manual", "all-pdf",
                             "LATEXMKOPTS=-interaction=nonstopmode -halt-on-error"])
            folder, basename = ROOT/"build/user-manual", "qdjj-user-reference"
            identity, suffix = "User's Guide and Reference Manual", "users-guide-reference"
            required = ("Installation and a first checked calculation", "Saved records and artifacts",
                        "Calculation input format", "Supported supplementary helpers")
        found, count = problems(folder, basename)
        pdf = folder/f"{basename}.pdf"
        found += artifact_problems(pdf, title=identity, version=version, pages=count, required=required)
        if found:
            raise RuntimeError(f"{manual} PDF validation failed:\n"+"\n".join(found))
        filename = f"qdjj-solver-{version}-{suffix}.pdf"
        shutil.copy2(pdf, output/filename)
        artifacts[manual] = dict(file=filename, sha256=digest(output/filename), pages=count)
        print(f"validated {filename}: {count} pages", flush=True)
    tools = {}
    for name in ("sphinx", "myst-parser", "sphinxcontrib-svg2pdfconverter", "pypdf", "matplotlib", "numpy", "scipy"):
        tools[name] = metadata.version(name)
    for name in ("xelatex", "pdflatex", "latexmk", "rsvg-convert"):
        tools[name] = subprocess.check_output([name, "--version"], text=True).strip().splitlines()[0]
    manifest = dict(format="qdjj-manuals", format_version=1, software_version=version, edition=label,
                    source_commit=commit, source_modified=dirty, tools=tools, manuals=artifacts)
    if args.source_archive:
        manifest["source_archive"] = dict(file=args.source_archive.name, sha256=digest(args.source_archive))
    (output/f"qdjj-solver-{version}-manuals.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
