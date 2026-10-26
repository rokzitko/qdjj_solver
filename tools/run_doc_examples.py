"""Execute canonical Markdown Python sessions in fresh processes and temporary directories.

python -B tools/run_doc_examples.py          # handbook chapters, one process per chapter
python -B tools/run_doc_examples.py --page docs/spectroscopy.md
"""

import argparse
import importlib.util
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[1]
THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
DIRECTIVE = re.compile(r"<!--\s*docs-test:\s*(continue|setup|skip|weekly|requires)\b(.*?)-->", re.S)


def code_blocks(relative, root=ROOT):
    steps, pending = [], None
    for token in MarkdownIt("commonmark").parse((root/relative).read_text()):
        if token.type in ("html_block", "html_inline"):
            for kind, argument in DIRECTIVE.findall(token.content):
                if kind in ("setup", "continue", "requires"):
                    steps.append((kind, argument.strip(), token.map[0]+1))
                else:
                    pending = kind
        elif token.type == "fence" and token.info.strip() == "python":
            steps.append((pending or "run", token.content, token.map[0]+1))
            pending = None
    return steps


def run_page(relative, namespace, executed, root=ROOT, require_all=False):
    if relative in executed:
        return
    executed.add(relative)
    have_tenpy = importlib.util.find_spec("tenpy") is not None
    weekly = os.environ.get("QDJJ_DOCS_WEEKLY") == "1"
    dependency = None
    for kind, code, line in code_blocks(relative, root):
        if kind == "requires":
            if code != "dmrg":
                raise ValueError(f"{relative}:{line}: unknown example dependency {code}")
            dependency = code
            continue
        if kind == "continue":
            run_page(code, namespace, executed, root, require_all)
            continue
        needed, dependency = dependency, None
        if kind == "skip" or (kind == "weekly" and not weekly):
            continue
        if not have_tenpy and (needed == "dmrg" or "dmrg" in code or "mps" in code):
            if require_all:
                raise RuntimeError(f"{relative}:{line}: DMRG dependency required for this example")
            continue
        try:
            exec(compile(code, f"{relative}:{line}", "exec"), namespace)
        except Exception as error:
            raise AssertionError(f"{relative}:{line}: {type(error).__name__}: {error}") from error


def handbook_pages(root=ROOT):
    text = (root/"docs"/"handbook"/"index.md").read_text()
    return [name+".md" for name in re.findall(r"^/(docs/[^\s]+)$", text, re.M)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page", help="one canonical Markdown session")
    parser.add_argument("--require-all", action="store_true", help="fail instead of skipping a missing backend")
    args = parser.parse_args()
    for variable in THREADS:
        os.environ[variable] = "1"
    if args.page:
        with tempfile.TemporaryDirectory(prefix="qdjj-doc-example-") as directory:
            os.chdir(directory)
            run_page(args.page, {"__name__": "__docs__"}, set(), require_all=args.require_all)
        print(f"passed {args.page}", flush=True)
    else:
        for page in handbook_pages():
            if not any(kind == "run" for kind, _, _ in code_blocks(page)):
                continue
            command = [sys.executable, "-B", str(Path(__file__).resolve()), "--page", page]
            if args.require_all:
                command.append("--require-all")
            subprocess.run(command, check=True, timeout=300, cwd=ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
