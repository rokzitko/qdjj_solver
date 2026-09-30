# Bath accuracy and calculation-time comparisons

This directory supports the [bath and reservoir guide](../../bath_representations.md).
The measurements compare cosh-transformed Gauss–Legendre quadrature, positive
surrogate fitting, and Gaussian quadrature directly in normal-state energy.
The last is supplied through `DiscreteBath` and labelled `linear-gl` in the
tables; it is not a separate bath-generation option in solver input files.

## Files and conventions

| File or directory | Contents |
|---|---|
| `input.json` | Initial set of parameter points, fitting settings, solver settings, and numbers of timing repetitions. |
| `refinement-input.json` | Additional reference-convergence calculations motivated by the initial MPS residuals. |
| `chain-input.json` | Independent physical-electron-chain spectrum refinements, calculated by `tools/refine_bath_chain.py`. |
| `finite-controls-input.json` | Unrestricted QP references for timing comparisons of the earlier star/chain DMRG routines. |
| `edge-refinement-input.json` | Additional QP-star MPS refinement for the weak-coupling gap-edge parameter point. |
| `manifest.json` | Initial calculation's solver identity, software versions, computer details, commands, and list of completed calculations. |
| `baths.json` | Explicit normal energies, positive contact weights, physical scales, and fitting settings. |
| `measurements.json` | Energies and observables, individual timing samples, state diagnostics, calculation settings, and failures from the initial study. |
| `kernel.csv`, `quadratic.csv` | Independent frequency-grid and finite-band noninteracting comparisons. |
| `many_body.csv`, `layouts.csv` | Initial interacting calculations and comparisons of the earlier star/chain DMRG routines. |
| `failures.json` | Failed construction or solver tasks; dependent tasks cannot run without a successful bath. |
| `refinement/` | Separately archived additional convergence measurements and their own input/manifest. |
| `chain/` | Electron-chain reference calculations, with full finite-Hamiltonian residual checks. |
| `controls/` | Small exact finite-bath comparisons for the earlier star/chain routines. |
| `edge/` | Independent QP-star spectrum reference, complementing the chain coordinates. |
| `summary.json` | Reference selection, remaining convergence indicators, unavailable tasks, and accuracy-target assessments. |
| `interacting_comparison.csv` | Differences from explicitly labelled empirical interacting references. |
| `coordinate_comparison.csv` | Isolated/coupled coordinates and truncations versus an unrestricted same-bath calculation. |
| `layout_comparison.csv` | Earlier star/chain routines' energy and spin errors against matching finite-bath QP results. |
| `refinement.csv` | Combined studies of QP cutoff, MPS bond dimension and ordering, and bath size, including unconverged calculations. |
| `fit_sensitivity.csv` | Effects of fitting-window, weighting, mesh-density, seed, and start-count choices. |
| `time_to_accuracy.csv` | Fastest **tested even-size** configuration per family at each accuracy target. |
| `tables.md`, `*.svg` | Tables and figures generated from the stored measurements. |

Energies are in units of $\Delta=1$. `levels` always counts signed spinful normal
levels **per reservoir**, so the cosh-grid `pairs` argument is `levels/2`.
`gamma` means total hybridization for the single-dot studies. Double-dot
couplings are stated separately for its two contacts.

`signed_gap` is $E_{\mathrm{odd}}-E_{\mathrm{even}}$. The current is reported as
$I/(2e\Delta/\hbar)=\langle\partial_\phi H\rangle/\Delta$; with $\Delta=1$,
its numerical value equals the stored `phase_derivative`. The double-dot paper
uses $e\Delta/\hbar$; multiply the benchmark current by two before comparing
to that paper's plotted normalization. Knight shift is
$\kappa=1-2\langle S^z_{\mathrm{dot}}\rangle$ in the spin-up doublet.

QP convergence checks and residuals with a finite `cutoff` concern the projected
Hamiltonian $P_qHP_q$. MPS checks concern the full finite Hamiltonian.
Neither establishes bath convergence. Blank CSV cutoff fields represent `null`
(full finite mode space), while blank bond fields indicate a QP calculation.

## Reading precision and time

The quadratic and weak-coupling references are independent continuum formulas.
Interacting references are empirical finite-bath refinements. The report records
their bath step, bond step, residual, and cross-family difference separately.
These quantities are not added into a purported rigorous continuum error bar.
A same-size competing bath can be less accurate than the reference; its
cross-family difference is therefore reported separately from the reference's
own refinement indicators.

For `reference_resolved=true` at target $\tau$, the selected reference must have
both a bath and bond comparison, their observable changes must be below
$\tau/5$, and its physical residual must be below $\tau/10$. This is an explicit
empirical reporting rule, not a theorem relating residual to observable error.
The lower-bond run can be a deliberately coarse, flagged calculation; its status
is retained. The selected MPS reference itself must satisfy its declared
finite-problem convergence checks. When only the largest finite-bath ED result
is available as a reference, it does not meet this continuum-resolution criterion.

The initial study takes five bath-construction samples and three many-body
samples after one warm-up. Supplemental reference calculations take one measured
sample after a warm-up and use two independent initial MPS trials; these establish
numerical references and are not used to rank the initial finite-ED runtimes.
Medians and interquartile ranges are descriptive statistics from these samples,
not confidence intervals. Fresh-fit totals add the measured median construction
time to the median reused-bath calculation time. Reused-bath calculations include
model construction, eigenstate calculation, observables, residuals, and numerical
diagnostics. Loading Python modules, evaluating reference integrals, and writing
results to disk are excluded.

Peak memory is the largest measured resident memory of the running program,
including initialization and warm-up, and is reported in bytes. The QP storage
limit instead applies to an estimate made before constructing the basis/matrix. Failed or
unconverged cases remain in the evidence. An empty accuracy-target result means
no tested configuration reaches the target relative to that reference. If
`reference_resolved=false`, no corresponding continuum-accuracy conclusion
follows. An empty entry does not establish a limit on what the method can achieve
with further refinement.

## Reproduce

Run from the main source directory after installing with `[dmrg,plots]`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B tools/benchmark_baths.py --profile full --dmrg --output results/bath-full
python -B tools/benchmark_baths.py \
  --input docs/benchmarks/baths/refinement-input.json \
  --profile full --dmrg --output results/bath-refinement
python -B tools/refine_bath_chain.py --output results/bath-chain
python -B tools/benchmark_baths.py \
  --input docs/benchmarks/baths/finite-controls-input.json \
  --profile full --output results/bath-controls
python -B tools/benchmark_baths.py \
  --input docs/benchmarks/baths/edge-refinement-input.json \
  --profile full --dmrg --output results/bath-edge
python -B tools/plot_bath_benchmarks.py results/bath-full \
  --supplement results/bath-refinement --supplement results/bath-chain \
  --supplement results/bath-controls --supplement results/bath-edge
```

Run the numerical commands one at a time when measuring elapsed time. Repeating
a command continues its unfinished calculations if the settings and code match
the original run. Completed measurements are saved individually. Baths are
identified by their coefficients and construction settings, independently of
fitting time. If an additional convergence study uses a bath from the original
study, its coefficients must match exactly; the plotting script retains the
original construction timings for those shared baths.

To regenerate the figures and tables from the archived measurements:

```sh
python -B tools/plot_bath_benchmarks.py docs/benchmarks/baths \
  --supplement docs/benchmarks/baths/refinement --supplement docs/benchmarks/baths/chain \
  --supplement docs/benchmarks/baths/controls --supplement docs/benchmarks/baths/edge
```

Use a new directory under `results/` for further calculations. Copy results here
after reviewing them for inclusion in the report. The source archive includes
the comparison scripts and data; the precompiled packages contain the solver.

The saved numerical-library information includes library names, versions,
processor architecture, and thread counts, with computer-specific directory
names removed.
