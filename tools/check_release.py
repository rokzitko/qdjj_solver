"""Check release version agreement without importing the package or native core."""

import argparse
import ast
from pathlib import Path
import re
import tomllib


def check(root, tag=None):
    project = tomllib.loads((root/"pyproject.toml").read_text(encoding="utf-8"))["project"]
    version = project["version"]
    tree = ast.parse((root/"src/qdjj_solver/__init__.py").read_text(encoding="utf-8"))
    runtime = next(ast.literal_eval(node.value) for node in tree.body
                   if isinstance(node, ast.Assign)
                   and any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets))
    citation = re.search(r"^version:\s*(\S+)\s*$", (root/"CITATION.cff").read_text(encoding="utf-8"), re.MULTILINE)
    if runtime != version or citation is None or citation[1].strip("\"'") != version:
        raise ValueError("pyproject.toml, __version__, and CITATION.cff versions must agree")
    if tag is not None and tag != f"v{version}":
        raise ValueError(f"release tag must be v{version}, got {tag!r}")
    print(f"Release version: {version}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag")
    args = parser.parse_args()
    check(Path(__file__).resolve().parents[1], args.tag)


if __name__ == "__main__":
    main()
