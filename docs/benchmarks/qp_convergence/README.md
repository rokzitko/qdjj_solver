# Evidence for the QP-cutoff convergence guide

This directory supports [the convergence guide](../../qp_convergence.md).
The main archive contains 12 projected QP calculations and three unrestricted
QP-number DMRG references for the symmetric, zero-field doublet at
$U/\Delta=2$, $\Gamma/\Delta=0.4$, $D/\Delta=100$, $\phi=3\pi/5$.
The QP calculations cover every cutoff from one to six at both 24 and 32
signed normal levels per reservoir.

## Contents

| File | Contents |
|---|---|
| `input.json` | Explicit settings for reproducing the study with the current package, plus the empirical reference-resolution floors. |
| `manifest.json` | Source-repository revision context, extraction-code hashes, input, record identifiers, and evidence checksums. |
| `baths.json` | Exact normal energies, weights, scales, and metadata, keyed by their SHA-256 fingerprints. |
| `measurements.json` | Measured physical observables, derived charge/spin values, QP weights, sectors, original solver settings, diagnostics, numerical environments, and source-file checksums. |
| `current/` | Fresh current-solver checks with their own input, manifest, exact baths, implementation hashes, and measurements. |
| `summary.json` | Observable-specific bath/bond changes, empirical resolutions, cutoff errors, and comparisons with the current solver. |
| `convergence.csv`, `table.md` | Full-precision long-form results and a readable generated table. |
| `observables.svg`, `errors.svg`, `weights.svg` | Figures regenerated from the numerical evidence. |

`P_0`, `P_1`, and `P_2` count physical dot electrons. `qp_weights[n]` counts bath
QPs. `C_d_bath` is **derived** using $3(q_d-P_1)/4$ for the spin-up SU(2)
doublet; independent direct-operator checks are in
[`examples/qp_convergence.py`](../../../examples/qp_convergence.py) and
[`tests/test_qp_convergence.py`](../../../tests/test_qp_convergence.py).
The scalar correlation concerns the total spin of both leads, not the
contact-weighted spin density.

## Source and precision

The original calculations accompany Teodor Iličin and Rok Žitko's manuscript
*Quasiparticle-resolved variational theory of Andreev spin qubits*, specifically
the maximum-QP-number study in
`ASQ_variational/qp_solver/benchmarks/asq_reference/REPORT.md`. Original JSON
files supply full-precision expectation values; no curve digitization or
extraction from rounded manuscript tables is involved.

Each compact record names and hashes its source file. The source Git revision
in the manifest is context; file hashes identify the exact working-tree
contents imported. The QP importer verifies model parameters, sectors,
residuals, and agreement among the original timing repetitions, retains the
first repetition's expectations, and records the maximum repetition spread.
The current package's reproduction settings live in `input.json`; each row's
`options` preserves the settings actually used for that measurement.

Four distinct bath records are retained. Differences between the two versions
at each bath size are quadrature-library roundoff, checked at relative tolerance
$2\times10^{-14}$ and absolute tolerance $10^{-14}$. Fresh cutoff ladders share
exactly the same bath object. The normal-state weights are never renormalized.

The legacy DMRG references have converged sweep flags, sweep energies, discarded
weights, normalization checks, and variance estimates. Their `residual` is
explicitly `null`: the current direct residual was not measured in those runs.
The fresh reference workflow measures the full residual and retains the
independent seed trials and sweep histories. The legacy reference did not fix
eta; the current run selects eta $+1$, reproducing the same bound doublet.
Both use the dark-vacuum reduction and remove the QP-number cutoff within it.

For each observable, the plotted reference resolution is

```math
\max\bigl(10^{-8},\;4(\lvert\delta_{\mathrm{bath}}\rvert+
\lvert\delta_{\chi}\rvert)\bigr).
```

The bath step compares $(N,\chi)=(24,128)$ with $(32,128)$, and the bond step
compares $(24,128)$ with $(24,192)$. The selected reference is $(24,192)$.
The formula is an explicit empirical reporting convention for each observable,
not a rigorous error estimate. Energy changes are in units of $\Delta$ and the
other changes are dimensionless. Nominal differences below the reference
resolution remain in CSV/JSON; open plot markers place them at that resolution.

## Reproduction and export

The complete commands are in the [guide](../../qp_convergence.md#reproducing-the-reference-study).
To regenerate this directory's figures and tables in an ignored working area:

```sh
python -B tools/plot_qp_convergence.py docs/benchmarks/qp_convergence \
  --verification docs/benchmarks/qp_convergence/current \
  --output results/qp-convergence-plots
```

To repeat the fresh current-solver evidence:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B tools/benchmark_qp_convergence.py --profile check --dmrg \
  --output results/qp-convergence-current
```

After reviewing new evidence, export its compact records with:

```sh
python -B tools/benchmark_qp_convergence.py --output results/qp-convergence-current \
  --export docs/benchmarks/qp_convergence/current
```

The optional **one-time migration** used to create the main archive is:

```sh
python -B tools/benchmark_qp_convergence.py \
  --import-paper ../ASQ_variational/qp_solver/benchmarks/asq_reference \
  --output results/qp-convergence-import
```

Normal solver runs, tests, analysis, and figure regeneration read only files in
this repository. The importer is available to audit the extraction when the
original manuscript checkout is available. Its output must be a new archive
directory, so a migration cannot silently mix records from different sources.

The source distribution includes all these assets and the scripts; the wheel
contains the installed solver packages. Routine tests validate archive
fingerprints, physical probability identities, empirical-reference selection,
and regenerated tables using the shipped data, and solve only small finite
models. The full accuracy study is an explicit command.
