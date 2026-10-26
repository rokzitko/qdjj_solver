# Documentation policy and tooling

The documentation has three parts, each addressing a different question.
Where topics overlap, one document gives the authoritative account and the others link to it.

| Part | Answers | Owns | Source |
|---|---|---|---|
| Physics manual | What is computed and what does it mean? | Models, equations, methods and their error control, physical interpretation, worked physics with figures | `docs/manual/` (LaTeX, PDF) |
| User and developer manual | How do I run, check, save and extend the code? | Installation, workflows, Python API, command line and input files, results and persistence, convergence practice, troubleshooting, reference tables, developer guide | Canonical Markdown in `README.md`, `CONTRIBUTING.md`, `docs/`; Sphinx builds HTML and a separate user/reference PDF with selected chapters |
| Tutorials | Show me, and let me run it. | Runnable calculations with checked expected results | `examples/NN_topic/` |

The research applications (`applications/`) reproduce published results; they follow
their own conventions, described in [applications](../applications/README.md).

## Definition of done

A change of features or behaviour is complete when it updates the physics manual, the user
manual and a tutorial, or the pull request states why one of them is unaffected. Tables of
options, results, command-line arguments, observables, probes, conventions and error
messages are generated from the code: update the source, never the generated page.

## Generated pages

`python tools/gen_docs.py` writes, and `python tools/gen_docs.py --check` verifies:

| Output | Source |
|---|---|
| `docs/reference/options.md`, `docs/reference/results.md` | Field docstrings of the option and result dataclasses |
| `docs/reference/cli.md` | `build_parser()` of `qdjj_solver.cli` and `qdjj_solver.qp_solver.cli` |
| `docs/reference/observables.md` | `OBSERVABLE_DOCS` and `PROBE_DOCS` |
| `docs/reference/conventions.md`, `docs/manual/generated/conventions.tex` | `docs/conventions.toml` |
| `docs/reference/concordance.md`, `docs/manual/generated/concordance.tex` | `docs/concordance.toml` |
| `docs/reference/api.md` | The public `__all__` lists (Sphinx autodoc; needs the `dmrg` extra) |
| `docs/reference/models.md` | `docs/parameters.toml`, callable signatures and model-class attribute docstrings |
| `docs/reference/configuration.md`, `docs/reference/records.md` | `docs/formats.toml`, checked against input/record contracts |
| Handbook code regions and `docs/manual/generated/practical_code.tex` | Marked runnable example excerpts, `examples/quickstart.py`, Tutorial 01's canonical JSON input and handbook reference values |
| `docs/troubleshooting.md` | `docs/troubleshooting.toml`; quoted messages must occur in the code |
| Tutorial table of `examples/README.md`, tutorial list of `docs/index.md`, generated sections of each tutorial README | `examples/catalog.toml` and the tutorials' `reference/*.json` |
| `docs/manual/generated/pointers.tex`, `docs/manual/generated/version.tex` | User-manual headings, `examples/catalog.toml`, `pyproject.toml` |

The numbers quoted in the physics manual, `docs/manual/generated/quoted.tex`, are written
separately by `make quoted` (see the [physics manual](manual/README.md)).

Generated files are committed so that they render on GitHub. Their first line names the
generator.

## Writing Markdown

Every shared page must preserve the same content and mathematical notation on GitHub and on the site:

- Use only GitHub Markdown. MyST directives appear only in the site-only pages
  `docs/index.md`, `docs/handbook/index.md` and `docs/reference/api.md`.
- Follow the mathematics rules in `AGENTS.md`: `$...$` inline, ``$`...`$`` for spans with
  escaped punctuation, and ```` ```math ```` fences for display equations. A multi-line display
  equation uses an explicit environment such as `aligned` and contains no blank lines.
- Keep mathematics out of headings, and link files rather than directories.
- Link other pages relatively. Absolute `https://github.com/rokzitko/qdjj_solver/blob/main/...`
  links are reserved for `README.md`, which is also the package description on PyPI; the
  site turns them into internal links.
- Refer to the physics manual by section title and label, for example
  “Probing the driven junction” (`ch:probes`). The checker verifies that the label exists and
  that a quoted title matches the LaTeX heading.

## Checks

| Command | When | What it checks |
|---|---|---|
| `python tools/check_docs.py` | Every pull request | Links and anchors; headings; physics-manual labels and quoted titles; generated files; docstrings; tutorials against the catalogue and the chapters' `\tutorial` pointers, and shared quantities; quoted numbers; listings and `\code{}` names against the package; hard-coded cross-reference numbers |
| `python -B -m pytest -m docs` | Every pull request | The checks above, negative tests that each check catches rot, and the Python code blocks of the user manual |
| `python -B -m pytest -m tutorial` | Every pull request | The quick profile of every tutorial against its reference |
| `python -m sphinx -W --keep-going -c docs -b html . build/docs` | Every pull request | Strict site build, including the API reference |
| `python -B tools/build_manuals.py --manual user` | Every pull request | User/reference PDF from shared sources, TeX diagnostics, actual PDF/version/content and manifest |
| `python -B tools/run_doc_examples.py --require-all` | Release source build | Independent handbook Python sessions, one fresh process and directory per chapter |
| `make -C docs/manual pdf`, then `python tools/check_manual.py` | Weekly | Physics manual without warnings or overfull boxes |
| `python tools/run_manual_listings.py` | Weekly | Every Python listing of the physics manual, in document order |
| `python tools/check_tutorials.py --profile full` | Weekly | The full profile of every tutorial |
| `make -C docs/manual data`, then `python -B docs/manual/scripts/quoted.py --compare OLD` | Weekly | Recomputed figure data leave the quoted numbers unchanged within their printed precision; `OLD` is a copy of `docs/manual/generated/quoted.json` taken before the recomputation |

The continuous-integration (CI) workflow `.github/workflows/tests.yml`, which GitHub runs for
every pull request, has a `docs` job. It runs the per-pull-request checks with fixed versions
of Sphinx and its extensions (`tools/docs-constraints.txt`) and uploads HTML and the
user/reference PDF. The weekly `docs-extended` job runs the remaining checks with
the latest versions and builds both PDFs. Neither manual has an overall page limit; page estimates for individual
additions are planning aids, not build constraints. The checker reports the page count
for information.

### Code blocks in Markdown

The Python blocks of the user-manual pages run in order, and later blocks use the variables
of earlier ones. HTML comments, invisible on GitHub and on the site, steer them:

- `<!-- docs-test: continue README.md -->` first runs the blocks of another page;
- `<!-- docs-test: setup CODE -->` runs code that the reader does not see, such as a
  frequency grid;
- `<!-- docs-test: skip -->` skips the next block, and `<!-- docs-test: weekly -->` runs it
  only in the weekly job.
- `<!-- docs-test: requires dmrg -->` declares the next block's optional backend;
  release checks require all declared dependencies rather than silently skipping them.

`tools/run_doc_examples.py` provides the common parser and example runner. New self-contained
examples show all required imports and inputs; continuation-style fragments
identify their preceding session. Marked excerpts from the standalone handbook
scripts are generated and tested through those full scripts, including CLI
equivalence and fresh-directory save/load checks. Skipping an excerpt is not a
substitute for testing its canonical executable source.

The [PDF publication guide](user_manual.md) defines document selection, typography,
offline links, system dependencies and the release-artifact contract. The release
job builds both manuals from the exact sdist, not a second checkout, and records
source and PDF checksums. Documentation tooling does not enter runtime dependencies
or installed wheels.

## Tutorials

Each tutorial lives in `examples/NN_topic/` and is listed in `examples/catalog.toml`, which
names its physics-manual labels (the first one's chapter prints `\tutorial{NN_topic}`), its
user-manual pages, and the quantities it shares with the physics manual. A tutorial has:

- `README.md` with the sections Goal, Physics in brief, Physics manual, Run, Expected results,
  What to look at, Try next, Files and See also; the physics-manual pointer and the table of
  expected results are generated;
- `run.py` with `--profile quick|full` and `--output DIR`, writing `results.json`;
- `reference/quick.json` and `reference/full.json` with values and tolerances;
- `expected.png`, the figure of the full profile (below 150 kB), when the tutorial has one.

To add one: write the directory, add its catalogue entry, add `\tutorial{NN_topic}` to the
chapter, and run `python tools/gen_docs.py` and the checks.

## Building

```sh
python -m pip install -e '.[dmrg,dev,docs]'
python -m sphinx -W --keep-going -c docs -b html . build/docs   # open build/docs/index.html
cd docs/manual && make pdf                                       # build/manual.pdf; see its README
```
