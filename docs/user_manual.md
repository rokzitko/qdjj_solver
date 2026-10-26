# PDF manuals: reading, building and publication

The **User's Guide and Reference Manual** explains how to perform complete calculations,
then provides constructor, solver, API, JSON and result-record references. Its
canonical Markdown pages and source descriptions are shared with the HTML site.
The **physics manual** develops the physical interpretation and methods, including
“From a Hamiltonian to a checked calculation” (`ch:practical`) before periodic driving.

## Download a matching edition

The [GitHub releases](https://github.com/rokzitko/qdjj_solver/releases) publish
the following assets for releases built with the two-manual workflow:

- `qdjj-solver-<version>-users-guide-reference.pdf`
- `qdjj-solver-<version>-physics-manual.pdf`
- `qdjj-solver-<version>-manuals.json`, recording source identity, tool versions
  and checksums of both PDFs and the source archive.

Earlier releases may predate PDF publication. Choose a manual matching the
installed software; development builds identify their source commit and local
changes. PDFs are generated artifacts, not files committed to the source tree.
The source distribution contains the canonical chapters, example inputs, scripts
and compact reviewed reference records for offline reproduction.

## Build locally

From the source root, install the documentation environment in your chosen
virtual environment:

```sh
python -m pip install '.[dmrg,dev,docs,docs-pdf,plots]' -c tools/docs-constraints.txt
python -B tools/build_manuals.py --manual all
```

The build script verifies that the versions agree and the generated sources are up to date, builds
from saved data, checks the actual PDFs and TeX diagnostics, then writes versioned
files and their manifest under `build/manuals/`. It does not run physics scans.
Use `--manual user` or `--manual physics` to build one document. Build logs are in
`build/manual-logs/`; renderer intermediates are in `build/user-manual/` and
`docs/manual/build/` respectively.

System requirements are XeLaTeX, pdfLaTeX, latexmk, biber, the Sphinx/physics
manual's TeX packages, TeX Gyre and DejaVu fonts, and `rsvg-convert` for vector
SVG figures. On Debian/Ubuntu, the CI environment installs:

```sh
sudo apt-get install texlive-xetex texlive-latex-extra texlive-bibtex-extra \
  texlive-science texlive-fonts-recommended fonts-dejavu-core \
  lmodern cm-super biber latexmk librsvg2-bin
```

On macOS, a full MacTeX installation supplies TeX and fonts; `rsvg-convert` is
provided by librsvg. On Windows, use a TeX distribution with XeLaTeX/latexmk,
the named fonts, Make and librsvg on PATH, or build in WSL. PDF CI targets
Linux; local macOS builds are also exercised during development.

The direct handbook commands are:

```sh
python -m sphinx -W --keep-going -t handbook -c docs -b latex . build/user-manual
make -C build/user-manual all-pdf LATEXMKOPTS='-interaction=nonstopmode -halt-on-error'
python tools/check_manual.py --build-dir build/user-manual --basename qdjj-user-reference --title "User's Guide and Reference Manual"
```

Use the build script when preparing files for publication. It also checks the required
content and records the build environment. Neither manual has a page ceiling.

## Maintain one source for each fact

- Edit workflow prose in the canonical `docs/` pages, not generated TeX.
- Edit option/result/model attribute descriptions in their Python sources.
- Edit constructor explanations in `docs/parameters.toml`; signatures and
  defaults are read from code and parameter membership is checked.
- Edit format explanations in `docs/formats.toml`; envelope contracts and
  representative emitted records are tested against the implementation.
- Edit shared runnable examples and their marked excerpts, then regenerate
  displayed code with `python tools/gen_docs.py`.
- Edit Tutorial 01's shared `input/reference.json` when changing its physical
  parameters, then recompute and review every affected reference result.

The handbook's explicit contents tree is `docs/handbook/index.md`. It includes
selected convergence evidence, with links to full reports. Printed reference
tables use labelled entries when columns would become unreadable. Both forms
use the same parsed content, and the PDF includes an index of public names and
parameters. Inspect rendered pages at reading size after layout changes.

## Offline links and release verification

Links within the handbook resolve internally. Links to excluded documents and
source/data files point to the source revision of the edition; they are not
HTML-only `_downloads` paths or implicit PDF attachments. Filenames identify
the corresponding files in the matching source archive.

The [release workflow](releases.md#versioned-artifacts) builds both PDFs from
the exact source distribution used to build the wheels. It runs small examples against the installed
package, verifies version identity and PDF completeness, and publishes the
validated files without rebuilding them. Reproduction means rebuilding the
document and checking its content; checksums identify a particular artifact,
not a promise of byte-identical TeX output across platforms.
