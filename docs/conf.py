"""Sphinx configuration for the qdjj_solver user and developer manual.

The source directory is the repository root, so that every Markdown page keeps
its GitHub-relative links. Build from the repository root with

    python -m sphinx -W --keep-going -c docs -b html . build/docs

Shared pages are plain GitHub Markdown; MyST directives appear only in the
site-only pages (``docs/index.md`` and ``docs/reference/api.md``).
"""

import sys
import os
import subprocess
import tomllib
from pathlib import Path

DOCS = Path(__file__).resolve().parent
sys.path.insert(0, str(DOCS/"_ext"))

import qdjj_solver  # noqa: E402  (the installed package, including the native extension)

project = "qdjj_solver"
author = "The qdjj_solver developers"
copyright = "the qdjj_solver developers"
release = qdjj_solver.__version__
version = release
source_version = tomllib.loads((DOCS.parent/"pyproject.toml").read_text())["project"]["version"]
if release != source_version:
    raise RuntimeError(f"installed qdjj_solver {release} does not match documentation source {source_version}")
handbook = bool(globals().get("tags") and globals()["tags"].has("handbook"))

extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.mathjax", "sphinx_copybutton",
              "github_markdown", "handbook"]

root_doc = "docs/index"
source_suffix = {".md": "markdown"}
# Sphinx translates "**/" into "one or more directories", so list both levels explicitly.
include_patterns = [
    "README.md", "CONTRIBUTING.md",
    "docs/*.md", "docs/**/*.md",
    "examples/README.md", "examples/*/README.md",
    "applications/README.md", "applications/*/README.md",
    "NRG_comparisons/README.md", "NRG_comparisons/test*/README*.md",
    "NRG_comparisons/test*/manual_run/README.md", "NRG_comparisons/test*/output/report.md",
    "NRG_comparisons/test4/FINAL_REPORT.md",
    "large_Gamma/README.md", "large_Gamma/output/report.md",
]
# Only exclusions prune the directory walk.
exclude_patterns = [
    ".git", ".venv", "build", "dist", "runs", "results", "reports", "**/runs", "**/__pycache__",
    ".hypothesis", ".ruff_cache", ".pytest_cache", "src", "cpp", "tests", "tools",
    "docs/manual/build", "docs/manual/figures", "docs/_root",
]

myst_enable_extensions = ["dollarmath", "gfm_autolink", "strikethrough"]
myst_fence_as_directive = ["math"]
myst_heading_anchors = 6

autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_class_signature = "separated"
autodoc_preserve_defaults = True
python_maximum_signature_line_length = 80

# A separate publication tree selects the PDF without changing the website's root.
# Keep the literal source lists above: tools/check_docs.py reads them without importing Sphinx.
if handbook:
    extensions.append("sphinxcontrib.rsvgconverter")
else:
    exclude_patterns.append("docs/handbook")

handbook_source_ref = os.environ.get("QDJJ_DOCS_REF", "")
if not handbook_source_ref:
    try:
        handbook_source_ref = subprocess.check_output(
            ["git", "-C", str(DOCS.parent), "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        handbook_source_ref = f"v{release}"
handbook_build_label = os.environ.get("QDJJ_DOCS_LABEL", f"Development snapshot {handbook_source_ref[:12]}")
today = handbook_build_label if handbook else ""
latex_documents = [("docs/handbook/index", "qdjj-user-reference.tex",
                    r"qdjj\_solver: User's Guide and Reference Manual", author.replace("_", r"\_"), "manual", False)] if handbook else []
latex_engine = "xelatex"
latex_use_xindy = False
latex_show_pagerefs = True
latex_show_urls = "no"  # source filenames and explicit versioned links are retained in print
latex_elements = {
    "releasename": "Version",
    "papersize": "a4paper", "pointsize": "11pt", "classoptions": ",oneside,openany",
    "cmappkg": "",
    "fontpkg": (r"\setmainfont{texgyretermes}[Extension=.otf,UprightFont=*-regular,ItalicFont=*-italic,BoldFont=*-bold,BoldItalicFont=*-bolditalic]"
                r"\setsansfont{texgyreheros}[Extension=.otf,UprightFont=*-regular,ItalicFont=*-italic,BoldFont=*-bold,BoldItalicFont=*-bolditalic]"
                r"\setmonofont{DejaVuSansMono}[Extension=.ttf,UprightFont=*,ItalicFont=*-Oblique,BoldFont=*-Bold,BoldItalicFont=*-BoldOblique,Scale=0.85]"),
    "preamble": (r"\providecommand{\lt}{<}\providecommand{\gt}{>}\setlength{\emergencystretch}{3em}"
                 r"\setlength{\headheight}{24pt}\hypersetup{bookmarksdepth=3}"
                 # Preserve qualified names in the two-column index, but permit
                 # a line break between module/class components.
                 r"\ExplSyntaxOn\cs_set_protected:Npn\sphinxstyleindexextra#1{"
                 r"\seq_set_split:Nnn\l_tmpa_seq{.}{#1}"
                 r"\space(\emph{\seq_use:Nn\l_tmpa_seq{.\allowbreak}})}\ExplSyntaxOff"),
    "sphinxsetup": "verbatimwithframe=true,verbatimwrapslines=true",
}
if handbook:
    # TeX Live supplies these files on macOS; Debian's system DejaVu package is
    # registered by family name instead. Resolve the installed font, not a substitute.
    located = subprocess.run(["kpsewhich", "DejaVuSansMono.ttf"], capture_output=True, text=True)
    if located.returncode:
        latex_elements["fontpkg"] = latex_elements["fontpkg"].split(r"\setmonofont", 1)[0] + (
            r"\setmonofont{DejaVu Sans Mono}[Scale=0.85]")

html_theme = "pydata_sphinx_theme"
html_title = f"qdjj_solver {release}"
html_extra_path = ["_root"]
html_theme_options = {
    "navigation_with_keys": False,
    "show_toc_level": 2,
    "secondary_sidebar_items": ["page-toc"],
}
html_sidebars = {"**": ["sidebar-nav-bs"]}
html_show_sourcelink = False

# The weekly link check reports external problems only; links into this repository are
# checked offline by tools/check_docs.py, and DOI resolvers often refuse automated clients.
linkcheck_ignore = [r"https://github\.com/rokzitko/qdjj_solver/", r"https://doi\.org/"]
linkcheck_anchors_ignore_for_url = [r"https://github\.com/"]
linkcheck_timeout = 30
