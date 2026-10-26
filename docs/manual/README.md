# Physics manual

This directory contains the LaTeX source of the `qdjj_solver` physics manual. The
manual describes the package from a physical perspective: the Hamiltonian, the equations
that are solved, the approximations, and the physical interpretation. It builds up the
physics of driven junctions step by step through worked examples with figures.

The opening chapters assume basic BCS theory and second quantisation, but introduce
superconducting impurities from the beginning: “Andreev reflection and subgap states”
(`ch:background`), “The model” (`ch:model`), “Magnetic impurities, Kondo physics and
screening” (`ch:correlations`), and “Subgap spectroscopy and Andreev qubits”
(`ch:experiments`). Original schematics, analytic limits and interpretation exercises
connect the Hamiltonian to tunnelling, microwave measurements and qubit experiments.

It is one of three documents ([documentation policy](../documentation.md)): the user and
developer manual explains how to run and extend the code, and each worked example has a
runnable [tutorial](../../examples/README.md). The last chapter, “From the equations to the
code” (`ch:equations-to-code`), maps the equations onto the code.

## Building

```sh
cd docs/manual
make pdf       # quoted numbers and figures from data/, then latexmk; output build/manual.pdf
make quoted    # generated/quoted.tex from data/*.json (no solver, no matplotlib)
make figures   # render figures/*.pdf from data/*.json only
make data      # recompute data/*.json with the solver (tens of minutes)
make clean     # remove build/ and figures/
python ../../tools/check_manual.py          # no warnings or overfull boxes; reports page count
python ../../tools/run_manual_listings.py   # every Python listing, in document order
```

Requirements:

- a TeX distribution with `pdflatex`, `latexmk`, `biber`, KOMA-Script and `biblatex-phys`;
- `matplotlib` (the `plots` extra) for the figures;
- the installed package for `make data`;
- `pypdf` (the `docs-pdf` extra) for checking the built PDFs with `check_manual.py`.

To build and validate both versioned PDF editions, use
`python -B tools/build_manuals.py --manual all` from the source root; see the
[PDF publication guide](../user_manual.md) for the complete environment and outputs.

The manual has no overall page limit. Its length can grow with the physics and examples
it needs to explain.

## Layout

| Path | Content |
|---|---|
| `manual.tex`, `preamble.tex` | Main document, packages, notation and boxed environments |
| `chapters/` | One file per chapter, plus the appendix |
| `references.bib` | Verified bibliography; alternative keys are kept in `ids` |
| `scripts/fig_*.py` | Figure scripts with `compute` (calculation → `data/`) and `plot` (`data/` → `figures/`) subcommands |
| `scripts/manual_common.py` | Shared models, provenance records and figure style |
| `scripts/quoted.py` | Collects `quoted(record)` of every figure script into `generated/quoted.tex` (`make quoted`) |
| `scripts/edition.py` | Writes transient build identity to `build/edition.tex`; release builds supply their tag |
| `data/` | Tracked numerical results, each recording the script, code version, git commit and environment |
| `generated/` | Committed generated input: quoted numbers, conventions and concordance tables, user-manual and tutorial pointers, version |
| `figures/`, `build/` | Rendered figures and LaTeX output; not tracked by git |

## Maintenance

- Update the manual together with every new feature: equations, interface, an example
  and, where needed, a figure script with tracked data.
- Plotting never runs the solver. Recompute the data explicitly with `make data` (or
  `python -B scripts/fig_<name>.py compute`) when the numerical behaviour changes.
- `scripts/fig_primer.py` computes analytic teaching examples, records their assumptions
  in `data/primer.json`, and renders the introductory plots and original vector schematics.
  Its current plot reuses `data/static.json`. Schematics are labelled as such and do not
  represent measured data; analytical curves and finite-bath solver results retain their
  distinct approximation scopes. Use `python -B scripts/fig_primer.py compute` to regenerate
  just these inexpensive analytic records, then `make pdf` to render them.
- Never type a number that comes from a calculation. Return it from the figure script's
  `quoted(record)` and write `\q{<figure>.<key>}` in the text; `make quoted` regenerates
  `generated/quoted.tex`, so a recomputation changes the text and the change is reviewed in
  the diff. A quoted qualitative claim (for example an inequality) can be guarded by raising
  an error in `quoted()` when the data no longer support it.
- Point to the user manual with `\usermanual{docs/page.md}` and to a tutorial with
  `\tutorial{NN_topic}`; both are generated, and unknown keys are LaTeX errors. Refer to
  figures, sections and equations with `\cref`, never with hard-coded numbers.
- Code that a listing needs but the text does not show (for example a frequency grid) goes
  into `%! code` lines directly before `\begin{lstlisting}`; a `%! skip` line excludes the
  listing. `tools/run_manual_listings.py` runs the listings in document order, each building
  on the earlier ones.
- `tools/run_manual_listings.py` also follows generated LaTeX inputs in document
  order. The practical chapter's setup is generated from Tutorial 01's canonical
  input, so its displayed parameters and executable listing remain synchronized.
- Cite original papers for methods and related applications. Add entries to
  `references.bib` only after checking their metadata against the publisher or Crossref.
