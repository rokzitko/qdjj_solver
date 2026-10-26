"""Documentation consistency: the repository passes tools/check_docs.py, and each check catches rot."""

import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("markdown_it")
pytestmark = pytest.mark.docs

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_docs", ROOT/"tools"/"check_docs.py")
check_docs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_docs)


@pytest.mark.parametrize("name", sorted(check_docs.CHECKS))
def test_repository_passes(name):
    assert check_docs.CHECKS[name](ROOT) == []


def tree(tmp_path, files):
    conf = ('include_patterns = ["README.md", "docs/*.md"]\n'
            'exclude_patterns = ["build"]\n')
    for name, text in {"docs/conf.py": conf, **files}.items():
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return tmp_path


def test_site_pages_follow_the_sphinx_patterns(tmp_path):
    root = tree(tmp_path, {"README.md": "# A\n", "docs/b.md": "# B\n", "docs/sub/c.md": "# C\n",
                           "notes.md": "# D\n"})
    assert check_docs.site_pages(root) == ["README.md", "docs/b.md"]


def test_links_and_anchors(tmp_path):
    root = tree(tmp_path, {
        "README.md": ("# Title\n\n## Getting started\n\n[ok](docs/b.md#second-part) [self](#getting-started)\n"
                      "[code](tools.py) [web](https://example.org/x)\n"
                      "[gh](https://github.com/rokzitko/qdjj_solver/blob/main/docs/b.md#second-part)\n"),
        "docs/b.md": "# B\n\n## Second part\n",
        "tools.py": "",
    })
    assert check_docs.check_links(root) == []
    (root/"README.md").write_text("# T\n\n[a](docs/b.md#third) [b](missing.md) [c](docs)\n[d](notes.md)\n"
                                  "[e](https://github.com/rokzitko/qdjj_solver/blob/main/docs/b.md#nope)\n"
                                  "![f](figure.svg)\n")
    (root/"notes.md").write_text("# Notes\n")
    problems = "\n".join(check_docs.check_links(root))
    for expected in ("missing anchor #third", "missing link target: missing.md", "directory links",
                     "not part of the site: notes.md", "missing anchor #nope", "missing image target"):
        assert expected in problems


def test_duplicate_headings_get_numbered_anchors(tmp_path):
    root = tree(tmp_path, {"README.md": "# T\n\n## Run\n\n## Run\n\n[x](#run-1) [`code` title](#code-title)\n\n"
                                        "## `code` title\n"})
    assert check_docs.check_links(root) == []


def test_headings_without_math(tmp_path):
    root = tree(tmp_path, {"README.md": "# Energy $E$\n"})
    assert check_docs.check_headings(root) == ["README.md:1: mathematics in a heading: Energy $E$"]


def test_physics_manual_labels(tmp_path):
    root = tree(tmp_path, {"README.md": "# T\n\nSee `ch:probes` and `sec:gone`; `x:y` is not a label.\n",
                           "docs/manual/chapters/a.tex": "\\chapter{P}\\label{ch:probes}\n"})
    assert check_docs.check_labels(root) == ["README.md:3: unknown physics-manual label `sec:gone`"]


def test_quoted_physics_titles_must_match(tmp_path):
    chapter = ("\\chapter{Probing the driven junction}\\label{ch:probes}\n"
               "\\section{Mollow triplet and Autler--Townes doublet}\\label{sec:mollow}\n"
               "\\section{The $g$-factor}\\label{sec:g}\n")
    root = tree(tmp_path, {
        "README.md": ("# T\n\nRead “Probing the driven\njunction” (`ch:probes`), “Mollow triplet and "
                      "Autler–Townes doublet” (`sec:mollow`)\nand “Anything” (`sec:g`).\n"),
        "docs/manual/chapters/a.tex": chapter})
    assert check_docs.check_labels(root) == []
    (root/"README.md").write_text("# T\n\nSee “Probing a junction” (`ch:probes`).\n")
    assert check_docs.check_labels(root) == [
        "README.md:3: `ch:probes` is titled “Probing the driven junction”, not “Probing a junction”"]


def test_command_line(capsys):
    assert check_docs.main(["links"]) == 0
    assert "passed" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        check_docs.main(["unknown-check"])


def copy_sources(tmp_path):
    """A copy of the sources that tools/gen_docs.py reads (no data, builds or research studies)."""
    import shutil
    skip = shutil.ignore_patterns(".git", "build", "dist", "runs", "results", "reports", "__pycache__", "*.so",
                                  "*.pyd", "applications", "large_Gamma", "NRG_comparisons", "benchmarks",
                                  "figures", "data", "tests", ".*cache", "*.png", "*.pdf")
    for name in ("src", "tools", "docs", "examples", "cpp"):
        shutil.copytree(ROOT/name, tmp_path/name, ignore=skip)
    for name in ("README.md", "CONTRIBUTING.md", "pyproject.toml"):
        shutil.copy(ROOT/name, tmp_path/name)
    return tmp_path


def test_generated_files_must_be_current(tmp_path):
    root = copy_sources(tmp_path)
    assert check_docs.check_generated(root) == []
    conventions = root/"docs"/"conventions.toml"
    conventions.write_text(conventions.read_text().replace("Units", "Unit system"))
    assert check_docs.check_generated(root) == [
        "docs/reference/conventions.md: stale generated file; run python tools/gen_docs.py",
        "docs/manual/generated/conventions.tex: stale generated file; run python tools/gen_docs.py"]


def test_field_docstrings_are_required_and_plain(tmp_path):
    root = copy_sources(tmp_path)
    source = root/"src"/"qdjj_solver"/"qp_solver"/"solver.py"
    text = source.read_text()
    text = text.replace('"""Number of lowest eigenstates computed in the selected sector."""',
                        '"""Number of states, see :func:`solve`."""')
    text = text.replace('    """Seed of the random starting vector of the iterative eigensolvers."""\n', "")
    source.write_text(text)
    problems = "\n".join(check_docs.check_docstrings(root))
    assert "SolverOptions.eigenpairs: use plain text" in problems
    assert "SolverOptions.seed has no field docstring" in problems


def test_observable_descriptions_cover_the_built_in_observables():
    from qdjj_solver import cosh_grid, reference_model
    from qdjj_solver.common.models import OBSERVABLE_DOCS
    junction = reference_model(cosh_grid(1, bandwidth=3.), symmetry=False, compress=False)
    assert set(junction.observables) == set(OBSERVABLE_DOCS)


def tutorial_tree(tmp_path, readme_sections=check_docs.TUTORIAL_HEADINGS):
    import json
    root = tree(tmp_path, {
        "examples/catalog.toml": ('[[tutorial]]\ndir = "01_demo"\ntitle = "Demo"\nsummary = "s"\n'
                                  'physics = ["ch:demo"]\nuser = ["README.md"]\nquick_seconds = 1\n'
                                  '[[tutorial.shared]]\nreference = "x"\nmanual = "demo.x"\nrtol = 1e-6\n'),
        "examples/01_demo/README.md": "# Tutorial\n\n" + "".join(f"## {s}\n\ntext\n\n" for s in readme_sections),
        "examples/01_demo/run.py": "",
        "examples/01_demo/reference/quick.json": json.dumps(dict(values=dict(x=dict(value=1.)))),
        "examples/01_demo/reference/full.json": json.dumps(dict(values=dict(x=dict(value=1.)))),
        "docs/manual/chapters/a.tex": "\\chapter{Demo}\\label{ch:demo}\n\\tutorial{01_demo}\n",
        "docs/manual/generated/quoted.json": json.dumps({"demo.x": dict(value=1., text="1", format="g")}),
        "README.md": "# Root\n"})
    return root


def test_tutorials_follow_the_catalogue(tmp_path):
    root = tutorial_tree(tmp_path)
    assert check_docs.check_tutorials(root) == []
    (root/"docs/manual/generated/quoted.json").write_text('{"demo.x": {"value": 1.1, "text": "1.1", "format": "g"}}')
    (root/"docs/manual/chapters/a.tex").write_text("\\chapter{Demo}\\label{ch:demo}\n\\tutorial{02_other}\n")
    (root/"examples/02_extra").mkdir()
    problems = "\n".join(check_docs.check_tutorials(root))
    for expected in ("catalog.toml lists", "add \\tutorial{01_demo}", "disagrees with the physics manual",
                     "\\tutorial{02_other} is not in the catalogue"):
        assert expected in problems


def test_tutorial_readme_sections(tmp_path):
    root = tutorial_tree(tmp_path, readme_sections=("Goal", "Run"))
    assert any("sections must include" in p for p in check_docs.check_tutorials(root))


# --- executed code blocks -----------------------------------------------------------------
# Python blocks of the user-manual pages run in order in one namespace (in a temporary
# directory). Hidden HTML comments steer them: `<!-- docs-test: continue PAGE -->` first runs
# the blocks of another page, `<!-- docs-test: setup CODE -->` runs hidden code at its position,
# `<!-- docs-test: skip -->` skips the next block and `<!-- docs-test: weekly -->` runs it only
# in the weekly job (QDJJ_DOCS_WEEKLY=1). DMRG blocks are skipped without TeNPy.
example_spec = importlib.util.spec_from_file_location("run_doc_examples", ROOT/"tools"/"run_doc_examples.py")
doc_examples = importlib.util.module_from_spec(example_spec)
example_spec.loader.exec_module(doc_examples)


def code_blocks(relative):
    """[(kind, code or argument, line)] for a page: kind is run, setup, skip, weekly or continue."""
    return doc_examples.code_blocks(relative)


def snippet_pages():
    pages = [p for p in check_docs.site_pages(ROOT)
             if p in ("README.md", "CONTRIBUTING.md") or p.startswith(("docs/", "examples/"))]
    return [p for p in pages if any(kind == "run" or kind == "weekly" for kind, _, _ in code_blocks(p))]


def run_page(relative, namespace, executed):
    doc_examples.run_page(relative, namespace, executed)


@pytest.mark.timeout(900)
@pytest.mark.parametrize("relative", snippet_pages())
def test_code_blocks_run(relative, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run_page(relative, {"__name__": "__docs__"}, set())


def test_code_block_directives(tmp_path, monkeypatch):
    text = ("# T\n\n<!-- docs-test: setup x = 2 -->\n\n```python\ny = x + 1\n```\n\n<!-- docs-test: skip -->\n\n"
            "```python\nraise RuntimeError\n```\n\n<!-- docs-test: weekly -->\n\n```python\nraise RuntimeError\n```\n")
    path = tmp_path/"example.md"
    path.write_text(text)
    assert [kind for kind, _, _ in doc_examples.code_blocks("example.md", tmp_path)] == ["setup", "run", "skip", "weekly"]
    namespace = {}
    monkeypatch.delenv("QDJJ_DOCS_WEEKLY", raising=False)
    doc_examples.run_page("example.md", namespace, set(), tmp_path)
    assert namespace["y"] == 3


def test_latex_code_and_references(tmp_path):
    import shutil
    root = tree(tmp_path, {
        "src/qdjj_solver/module.py": "def solve_floquet(periodic, harmonics=4):\n    return periodic.quasienergies\n",
        "docs/manual/manual.tex": "\\input{chapters/a}\n",
        "docs/manual/chapters/a.tex": ("As in Fig. 9.1, \\code{solve_floquet} and \\code{harmonics=4} work, "
                                       "\\code{solve_magic} does not; see \\code{docs/manual/manual.tex}.\n"
                                       "\\begin{lstlisting}\nfrom qdjj_solver import NotThere\n"
                                       "states = solve_floquet(p, harmonics=4, magic=1)\nstates.quasienergies\n"
                                       "states.vanished\n\\end{lstlisting}\n")})
    (root/"tools").mkdir(exist_ok=True)
    shutil.copy(ROOT/"tools"/"run_manual_listings.py", root/"tools"/"run_manual_listings.py")
    problems = "\n".join(check_docs.check_latex(root))
    for expected in ("hard-coded reference 'Fig. 9.1'", "\\code{solve_magic} names nothing", "has no NotThere",
                     "keyword magic= does not occur", "attribute .vanished does not occur"):
        assert expected in problems
    assert "solve_floquet} names" not in problems and "manual.tex} is not a path" not in problems
