# Probe matrix elements and spectra

Calculating a spectrum requires the probe operator and initial-state populations, as well as the energies.
Specify which operator acts, which sector contains the final states and which
initial states are occupied. The physics manual develops the selection rules in
“Subgap spectroscopy and Andreev qubits” (`ch:experiments`) and
“Probing the driven junction” (`ch:probes`).

## Start with a static matrix element

The [static workflow](static_calculations.md) computes both the pair matrix
element and electron-addition strength. Rows are final bra states and columns
are initial ket states. A charge/phase probe preserves parity; an electron
operator changes the probed junction's parity. The complete system, including
the electron source, still conserves parity.

`matrix_elements` returns amplitudes. Combine coherent terms with their complex
phases intact, then square the resulting magnitudes to obtain strengths.
Energy differences determine possible line positions, while matrix elements
and populations determine visibility.

## A small independent Floquet spectrum

The following complete example drives a proximitized dot. It is deliberately
small so that all final states in its even sector can be retained.

```python
import numpy as np
from qdjj_solver import Hamiltonian, PeriodicHamiltonian, Sector
from qdjj_solver import create, hermitian_pair, number
from qdjj_solver.qp_solver import solve_floquet

def paired_dot(detuning=.2):
    n = number(0) + number(1)
    op = detuning*(n-1) - .3*hermitian_pair(create(0)*create(1))
    return Hamiltonian(op, 2, (1, -1), observables={"impurity_charge": n})

periodic = PeriodicHamiltonian.from_callable(
    lambda t: paired_dot(.2 + .02*np.cos(.72*t)), omega=.72,
    reference=paired_dot(), assume_fixed_coordinates=True)
states = solve_floquet(periodic, cutoff=0, sector=Sector(0, 0),
                       options=dict(floquet_cutoff=4, selection="all", method="dense"))
table = states.matrix_elements("charge")
initial = states.dressed_state(0)
lines = table.lines(initial=initial)
grid = np.linspace(-1.5, 1.5, 1001)
spectrum = table.spectrum(grid, initial=initial, width=.01, kind="signed")
print(lines.frequencies, lines.strengths)
assert np.all(np.isfinite(spectrum))
```

The explicit assumption of fixed coordinates is appropriate here because the specified
modes are fixed physical dot modes at every detuning. Use this option after
verifying that the canonical coordinates remain fixed throughout the drive.

## Inspect the line list before broadening

Floquet transition frequencies include integer multiples of the drive frequency.
`FloquetTransitions.photons` holds those harmonic indices;
`elements`, `frequencies` and `strengths` share the documented bra/ket axes.
`lines` merges degenerate multiplets, making the combined strength invariant under
changes of basis within a retained degenerate subspace. It omits the same-state,
zero-photon expectation contribution while retaining elastic harmonics.

The named probes and their physical meanings are listed in the
[observable/probe reference](reference/observables.md). Under a parameter drive,
a named observable uses its sampled time-dependent harmonics. Replacing a driven
current operator by its static value changes the probe problem.

## Choose populations and line-shape conventions

- `kind="signed"` displays absorption at positive frequency and stimulated
  emission at negative frequency.
- `kind="rate"` combines the positive- and negative-frequency contributions at
  a positive probe frequency.
- `kind="net"` subtracts stimulated emission from absorption.
- `width` is the half width at half maximum of the selected line shape.

Supply the initial-state populations explicitly through their weights. Predicting
steady-state populations requires an additional model of relaxation. Broadening
specifies the plotted line shape; physical line widths require a dynamical model
of the relevant decay and dephasing processes. A visible line requires both a
nonzero matrix element and population in its initial state.

## Check captured strength and convergence

Inspect `completeness`, the fraction of the sum-rule strength included in the calculation.
Increase the number of final states and the harmonic window, and verify that the spectral lines are stable. A complete
finite basis still needs QP and bath convergence. For transitions between different
sectors, use results from the same periodic Hamiltonian and scalar gauge.

For a realistic junction, first check the static spectrum and weak-drive limit
using a small retained state space. Then use [Tutorial 08](../examples/08_two_tone/README.md),
[photon-assisted tunnelling](../examples/09_photon_assisted_tunnelling/README.md)
and [spin resonance](../examples/10_spin_resonance/README.md) for worked extensions.
The [Floquet guide](floquet.md#transitions-and-spectra) describes compatibility
rules and the [result reference](reference/results.md) gives complete array layouts.
