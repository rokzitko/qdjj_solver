# Tutorial 00: a first calculation

## Goal

Reproduce the [main README's first calculation](../../README.md#first-calculation): the lowest
odd-parity state of a small superconducting junction with the QP expansion. Then check it
against DMRG and exact diagonalization of the same finite Hamiltonian.

## Physics in brief

A quantum dot with charging energy $U$ is coupled to two superconducting leads with a phase
difference $\phi$. Each lead is represented by a finite set of BCS levels (here one level pair
per lead). In the odd-parity sector the lowest state is the spin doublet $D$, in which the dot
retains a local moment that the quasiparticles of the leads partly screen. The QP expansion
limits the number of bath quasiparticles (the QP cutoff). Two quasiparticles already give
the energy to within about $10^{-3}$ of the exact finite-bath value. The probability of
finding no bath quasiparticle measures how strongly the leads dress the dot.

## Physics manual

<!-- BEGIN generated: physics-manual (tools/gen_docs.py) -->
The physics behind this tutorial is explained in the
[physics manual](../../docs/manual/README.md):

- section “A first calculation” (`sec:first-calculation`)
- section “The quasiparticle expansion” (`sec:qp`)
<!-- END generated: physics-manual -->

## Run

From the repository root, after `python -m pip install .`:

```sh
python -B examples/00_first_calculation/run.py                   # quick: QP only, < 1 s
python -m pip install '.[dmrg]'
python -B examples/00_first_calculation/run.py --profile full    # adds DMRG and exact, a few seconds
```

The quick profile performs the QP calculation and checks the energy, residual and
normalization asserted in the main README. The full profile also checks that the three
calculations use the same Hamiltonian and sector and that DMRG agrees with the exact energies.
Results are saved under `runs/examples/00_first_calculation/<profile>/`: `qp.json`, and in the
full profile `dmrg.json`, `exact.json` and `comparison.json`, together with the summary
`results.json`. Use `--output DIRECTORY` for another destination and `--save-states` to save
the wavefunctions.

Each input can also be solved directly from the command line (the DMRG input needs the
`dmrg` extra):

```sh
qdjj-solver solve examples/00_first_calculation/input/qp.json --output runs/first/qp.json
qdjj-solver solve examples/00_first_calculation/input/dmrg.json --output runs/first/dmrg.json
qdjj-solver solve examples/00_first_calculation/input/exact.json --output runs/first/exact.json
```

## Expected results

<!-- BEGIN generated: expected-results (tools/gen_docs.py) -->
`python tools/check_tutorials.py 00_first_calculation` runs the quick profile and compares it with
[`reference/quick.json`](reference/quick.json); add `--profile full` for
[`reference/full.json`](reference/full.json).

| Quantity | Key | Quick profile | Full profile | Tolerance |
|---|---|---|---|---|
| Lowest odd-parity energy at QP cutoff two | `qp_energy` | -0.295059596719 | -0.295059596719 | ±1e-10 |
| Eigensolver residual at QP cutoff two | `qp_residual` | at most 1e-09 | at most 1e-09 | upper bound |
| Probability of no bath quasiparticle | `qp_weight_0` | 0.89270561 | 0.89270561 | ±1e-08 |
| Exact finite-bath energy of the lowest state | `exact_energy_0` | – | -0.296312491144 | ±1e-09 |
| Exact finite-bath energy of the second state | `exact_energy_1` | – | 1.727646371364 | ±1e-09 |
| Largest difference between the DMRG and exact energies | `dmrg_minus_exact` | – | at most 1e-09 | upper bound |
| Cutoff-two energy above the exact value | `qp_minus_exact` | – | 0.001252894425 | ±1e-09 |
<!-- END generated: expected-results -->

The archived result records in [`reference/`](reference/qp.json) use the standard
`qdjj-eigenstates` JSON format: energies, residuals, observables, numerical settings, the
complete Hamiltonian and bath records, software versions and checksums that identify the
solver code. [`reference/comparison.json`](reference/comparison.json) records the
DMRG-minus-exact energy differences, the cutoff-two shift of the lowest energy and the DMRG
matrix of the impurity spin printed by the main README. The rows and columns of this matrix
are the bra and ket states in order of energy, the spin is in units of $\hbar$, and complex
arrays are stored as `real` and `imag` fields.

## What to look at

- The QP weights, approximately `[0.8927, 0.1010, 0.0063]` for zero, one and two bath
  quasiparticles: the state is mostly a dot spin dressed by a few quasiparticles.
- The cutoff-two energy lies $1.25\times10^{-3}$ above the exact value; the DMRG and exact
  energies agree to rounding error.
- The impurity spin of the doublet, the first diagonal element of the DMRG spin matrix in
  `comparison.json`, is $0.444$: the lead quasiparticles carry the rest of the total
  $S_z=1/2$.
- The residuals measure eigensolver accuracy within the chosen space; convergence in the
  QP cutoff, bond dimension and bath discretization requires separate checks.
- The lowest state of the selected sector is not necessarily the global ground state; to find
  it, compare the even and odd sectors.

## Try next

`run.py` asserts the reference values, so make these changes in copies of the inputs and
solve them with `qdjj-solver solve`, as shown above.

- Raise the QP `cutoff` in `qp.json` from 2 to 3 and watch the energy approach the exact
  value.
- Change `phi` or `u` in all three inputs and follow the doublet energy.
- Solve the even sector (`"parity": 0, "twice_sz": 0`) and compare its lowest energy with the
  doublet to find the ground-state parity.

## Files

| Input | Calculation | Requested states | Archived output |
|---|---|---|---|
| [`input/qp.json`](input/qp.json) | QP cutoff two | 1 | [`reference/qp.json`](reference/qp.json) |
| [`input/dmrg.json`](input/dmrg.json) | DMRG, `chi_max=32`, convergence required | 2 | [`reference/dmrg.json`](reference/dmrg.json) |
| [`input/exact.json`](input/exact.json) | QP backend, `cutoff=null` (full finite space) | 2 | [`reference/exact.json`](reference/exact.json) |

All inputs specify the same symmetric two-lead junction of `reference_model`
(`"kind": "reference"`), without compression. Its parameters differ from those of the
reference junction of the physics manual:

- `cosh-grid`, `pairs=1`, gap $\Delta=1$, normal-state **half-bandwidth** $D=5$;
- charging energy $U=1.6$, total hybridization $\Gamma=0.3$ ($0.15$ per lead);
- phase difference $\phi=0.8$ radians, with lead phases $[-\phi/2,+\phi/2]$;
- default zero detuning and zero field;
- odd parity and `twice_sz=1`, that is $S_z=1/2$;
- `compress=false`, seed `1729`, and one thread.

Energies use the gap as the unit and the solver's common
[energy reference](../../docs/models.md#units-and-energy-reference). The QP input requests one
state to match the main README; the command-line example [`examples/qp.json`](../qp.json)
requests two. [`reference/quick.json`](reference/quick.json) and
[`reference/full.json`](reference/full.json) hold the values and tolerances checked by
`python tools/check_tutorials.py 00_first_calculation`.

To regenerate the archived records, from the repository root with the DMRG extra installed:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B examples/00_first_calculation/run.py --profile full --output examples/00_first_calculation/reference
```

This runs the numerical checks before saving; review the numerical changes and the recorded
provenance, and remove the `results.json` written next to the records.

## See also

- [Using the solver and saving results](../../docs/interfaces.md), including
  [saved results and reproducibility](../../docs/interfaces.md#saved-results-and-reproducibility)
- [Numerical methods and convergence](../../docs/numerics.md)
- The [tutorial index](../README.md)
