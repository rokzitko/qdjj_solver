# Research applications

Ten worked investigations of published superconducting Anderson-model results.
Each directory contains the physical parameters, literature references,
calculation and plotting scripts, numerical results, and a scientific report.

| Application | Published result | Physics and numerical methods |
|---|---|---|
| [Impurity Knight shift](zitko_2023_knight_shift/README.md) | Pavešić et al., SciPost Phys. **15**, 070 (2023), Figs. 9(b), 10(a) | Doublet Zeeman response, spin screening, phase dependence, perturbation theory and NRG |
| [Double-dot supercurrent](zonda_2023_double_dot/README.md) | Žonda et al., PRB **107**, 115407 (2023), Fig. 9(a) | Two interacting orbitals, parity transitions, Josephson current, GAL and NRG |
| [Multi-terminal hidden symmetry](zalom_2024_multiterminal/README.md) | Zalom et al., PRL **132**, 126505 (2024), Fig. 2(d) | Three reservoirs, exact geometric mapping, subgap degeneracy and current conservation |
| [Surrogate subgap spectrum](paaske_2023_surrogate_spectrum/README.md) | Baran, Frost and Paaske, PRB **108**, L220506 (2023), Fig. 3 | Published bath coefficients, even/odd bath-size effects, coupled excitations and screening |
| [Gate–phase parity diagram](bargerbos_2022_parity_diagram/README.md) | Bargerbos et al., PRX Quantum **3**, 030311 (2022), Fig. S1(c) | Gate tuning away from half filling, deposited NRG data, charge response |
| [Kondo/Josephson crossover](choi_2004_kondo_josephson/README.md) | Choi, Lee, Kang and Belzig, PRB **70**, 020502 (2004), Fig. 3(a–c) | Gap/Kondo-scale competition, $0\text{–}\pi$ transitions, historical NRG comparison |
| [Unequal superconducting gaps](zonda_2016_unequal_gaps/README.md) | Žonda et al., PRB **93**, 024523 (2016), Fig. 8 | Distinct reservoir gaps, gate-dependent current, lead-current conservation and diagrammatic/NRG comparison |
| [Contact asymmetry](kadlecova_2017_asymmetry/README.md) | Kadlecová et al., PRB **95**, 195114 (2017), Fig. 1 and Eqs. (5)–(7) | Gate–phase boundaries, exact asymmetry mapping and current Jacobian |
| [Near-gap spectral peaks](hecht_2008_gap_edge/README.md) | Hecht et al., JPCM **20**, 275213 (2008), Fig. 6 and Fig. 7 diagnostic | Quadratic bath resolvent, many-body Lehmann/Lanczos spectra, spectral sum rules and explicit continuum-resolution limits |
| [Padé chain expansion and GAL](bobok_2025_chain_expansion/README.md) | Bobok, Frk, Pokorný and Žonda, PRB **112**, 205418 (2025), Figs. 8, 17 | Finite-/wide-band ChE baths, GAL currents, coherent common-lead coupling, same-parity singlet–triplet transitions and spin correlations |

## Running a case

The separate [large-hybridization benchmark](../large_Gamma/README.md) compares
2QP, 4QP, 6QP and DMRG across eight interactions and two phases. It has its own
restartable scan and analyzer; it extends the solver-convergence studies rather
than reproducing a particular paper figure.

From the main directory of the downloaded source code, install the solver:

```sh
python -m pip install '.[dmrg,plots]'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -m applications.zitko_2023_knight_shift.run --profile quick
```

Replace the case name to run another application. A **profile** selects a
predefined set of parameter points and numerical settings. All
`quick` and `paper` profiles use the base NumPy/SciPy/QP installation. The double-
dot convergence study and the Choi case's separate DMRG comparison require
TeNPy, as does the Bobok case's physical-chain comparison. Plotting requires
installation with `[plots]`; generating numerical data does not require Matplotlib.

The calculation scripts accept the following common options:

```sh
python -B -m applications.<case>.run --profile quick
python -B -m applications.<case>.run --profile paper --output results/my-case
python -B -m applications.<case>.run --profile convergence
python -B -m applications.<case>.plot --output results/my-case
```

`<case>` means one of the directory names above. New results are saved in
`applications/<case>/runs/<profile>/`. The results discussed in the reports
are archived in `output/quick/`, `output/paper/`, and `output/convergence/`.

- **quick:** small finite-bath examples with reference values for numerical checks.
- **paper:** the stated published curves or phase map, plus independent checks.
- **convergence:** selected points for bath-size, fit-window, and many-body convergence studies;
  the spectral case also varies real-frequency resolution and broadening.

`--input` selects another application JSON input file, in format version 1.
These inputs describe a parameter scan. The calculation scripts construct the
models and call the same solver functions used in the smaller examples.
`--raw` additionally saves each solution in the standard `qdjj-eigenstates`
JSON format, including its complete Hamiltonian and the identity of the solver
code used.

## Reproducibility and interpretation

Energies use **the specified gap as unit** ($\Delta=1$, or $\Delta_L=1$ for unequal
gaps), and `bandwidth` is the normal-state **half-bandwidth**.
Lead weights represent the normal-state measure, not the divergent
BCS quasiparticle density of states. Surrogate reservoirs use the method of
[Baran, Frost and Paaske, PRB **108**, L220506 (2023)](https://doi.org/10.1103/PhysRevB.108.L220506).
They are fitted once per distinct specification and reused throughout each scan,
or reconstructed from the authors' deposited coefficients in the surrogate-spectrum
case. The Choi case converts published bandwidth-unit parameters to gap units
separately for each $\Delta/T_K$; each resulting bath is archived. The Bargerbos
case sweeps detuning explicitly and keeps a fixed, uncompressed impurity basis.
The unequal-gap case fits each reservoir at its own gap. The spectral case uses
positive logarithmic quadrature for its high-resolution quadratic calculation
and surrogate baths for its interacting finite-space calculation; it converts
back to the paper's $D=1$ units when presenting spectra.
The Bobok case also uses [Padé chain-expansion baths](../docs/chain_expansion.md).
Its common-lead calculation keeps the full cross-hybridization matrix and both
interacting dot orbitals, and resolves total spin including the reservoir.

Each output directory includes:

- CSV tables of physical observables and comparison/convergence quantities;
- `eigenstates.json`: sector-resolved energies, expectation values, and numerical diagnostics;
- `manifest.json`: a calculation record containing the complete input, explicit
  bath nodes/weights, software versions, computer details, checksums identifying
  the solver and calculation scripts, and elapsed time;
- `summary.json`: calculated error metrics and important results;
- `comparison.svg` for the paper profile, generated from the archived tables.

The plotting scripts also save a PNG image. After reproducing the archived
`output/` directories, `python -B -m applications.analyze` updates each case's
`output/validation.json` with the convergence differences quoted in its report.

To reuse a bath without refitting it, copy its coefficients and settings from
`manifest.json` into the relevant `bath_record` entry of the input file. The
small numerical reference checks do this automatically. `cutoff` omitted or
`null` means unrestricted finite-space QP diagonalization; a finite cutoff is
always stated explicitly.
Single-dot dark-mode compression is documented and independently checked.

The case READMEs distinguish finite-Hamiltonian solver accuracy, bath convergence
and agreement with literature. Published rounded values and graphically extracted
NRG points have less precision than an eigensolver residual. Optional extraction
scripts retrieve literature data; the stored reference tables suffice for
ordinary calculations and tests, which require no network access.
The Choi case reproduces the key competing phases but retains sizeable numerical
differences from the 2004 curves at the stated parameters; its report quantifies
these rather than adjusting the coupling or current normalization. Reference
extraction for the Paaske case uses optional PyMuPDF. The Bargerbos extractor
downloads only the needed numerical data file from a large ZIP archive.
The Hecht report reproduces the quadratic near-edge peaks and demonstrates why
the interacting continuum peaks remain unresolved at the studied finite bath
sizes. Its archived spectral moments and small ground-state residuals do not
establish continuum convergence.

## Numerical checks and further studies

```sh
python -B -m pytest tests/test_applications.py tests/test_applications_second_set.py tests/test_applications_third_set.py
python -B -m pytest tests/test_applications.py tests/test_applications_second_set.py tests/test_applications_third_set.py -m slow
python -B -m pytest tests/test_chain_expansion.py tests/test_application_chain_expansion.py
python -B -m pytest tests/test_application_chain_expansion.py -m slow
```

New studies should follow the same directory structure and identify a specific
published result, its independent reference, all parameter conversions, and a
small calculation suitable for a routine numerical check. Keep downloaded
archives and large saved wavefunctions in `runs/`. The source archive includes
these applications; install the solver normally and run the scripts from the
main source directory.
