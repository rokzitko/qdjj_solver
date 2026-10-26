# Build the physics manual.
#   make pdf      quoted values and figures from the tracked data, then the PDF (build/manual.pdf)
#   make quoted   numbers quoted in the text, generated/quoted.tex from data/*.json (no solver)
#   make figures  render figures/*.pdf from data/*.json (fast; matplotlib)
#   make data     recompute data/*.json with the solver (minutes; explicit)
#   make clean    remove build artifacts and rendered figures
PYTHON ?= python -B
FIGURES := static primer dot tls highfreq numerics phase spectrum gate esr twotone pat esr_probe

.PHONY: pdf edition quoted figures data clean
pdf: edition quoted figures
	latexmk manual.tex

edition:
	$(PYTHON) scripts/edition.py

quoted:
	cd scripts && $(PYTHON) quoted.py

figures:
	cd scripts && for name in $(FIGURES); do $(PYTHON) fig_$$name.py plot || exit 1; done

data:
	cd scripts && for name in $(FIGURES); do $(PYTHON) fig_$$name.py compute || exit 1; done

clean:
	latexmk -C manual.tex
	rm -rf build figures
