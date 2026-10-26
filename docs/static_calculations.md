# Static states, observables and transitions

The aim is to find the ground-state parity, the lowest parity-changing and pair
excitations, and the current of a superconducting dot junction. This chapter uses
the reference junction of the physics manual and Tutorial 01. Its finite basis is larger than
that of the installation example, but the calculation remains inexpensive at QP cutoff two.

## Declare the physical and numerical problem

The shared input is
[`examples/01_subgap_states/input/reference.json`](../examples/01_subgap_states/input/reference.json).
It is also a complete unified-CLI input for the even sector. Its physical model
and bath parameters are also used by the tutorial and physics-manual figure scripts.

<!-- BEGIN generated: reference-input (tools/gen_docs.py) -->
```json
{
  "format_version": 1,
  "backend": "qp",
  "bath": {
    "kind": "cosh-grid",
    "pairs": 4,
    "delta": 1.0,
    "bandwidth": 10.0
  },
  "model": {
    "kind": "reference",
    "u": 1.0,
    "gamma": 0.5,
    "phi": 0.8,
    "detuning": 0.2,
    "field": 0.0,
    "symmetry": false,
    "compress": false
  },
  "sector": {
    "parity": 0,
    "twice_sz": 0
  },
  "cutoff": 2,
  "solver": {
    "eigenpairs": 2
  }
}
```
<!-- END generated: reference-input -->

`model` and `bath` define a finite physical Hamiltonian. `sector`, `cutoff` and
`solver` define the calculation performed on it. In particular, changing the
cutoff changes the retained many-body space at fixed bath; bath refinement changes
the finite Hamiltonian approximating the continuum. The labels and units are explained in
[physical models and conventions](models.md).

## Run the complete example

Install the base package, then run from the source root:

```sh
python -B examples/handbook/static.py --output runs/handbook/static
qdjj-solver solve examples/01_subgap_states/input/reference.json --output runs/handbook/cli-even.json --save-states
```

The [complete Python script](../examples/handbook/static.py) reads that JSON,
constructs the model, solves even and both odd spin sectors, evaluates matrix
elements, and verifies that saving and reloading the states preserves the results. The CLI command computes the
matching even states; use the Python script for the full cross-sector analysis.

The relevant part of the tested script is shown below. `config` is the parsed input
above and `output` is the directory selected by `--output`; imports and argument
handling are included in the downloadable script.

<!-- BEGIN generated: static-example (tools/gen_docs.py) -->
<!-- docs-test: skip -->
```python
h = model_from_config(config, default_compress=False)
even = solve(h, cutoff=config["cutoff"], sector=Sector(0, 0), options={"eigenpairs": 2})
odd = solve(h, cutoff=config["cutoff"], sector=Sector(1, 1))
partner = solve(h, cutoff=config["cutoff"], sector=Sector(1, -1))
signed_gap = float(odd.energies[0] - even.energies[0])
pair_gap = float(even.energies[1] - even.energies[0])
addition = odd.matrix_elements("d_up_dagger", even)
pair = even.matrix_elements("current")
current = float(np.real(even.observables["phase_derivative"][0]))
save_result(even, output/"even.json", save_states=True)
restored = load_result(output/"even.json")
np.testing.assert_allclose(restored.energies, even.energies, rtol=0, atol=1e-13)
```
<!-- END generated: static-example -->

## Read the energies before plotting

<!-- BEGIN generated: static-results (tools/gen_docs.py) -->
| Quantity | Expected value | Meaning |
|---|---:|---|
| `signed_gap` | 0.09578810633 | Signed doublet-minus-singlet energy, in gap units |
| `pair_gap` | 0.7159930829 | Even pair excitation, in gap units |
| `current` | 0.05034773273 | Singlet current in units of 2e Delta/hbar |
| `ground_parity` | 0 | Global ground-state parity (zero is even) |
| `addition_weight` | 0.4821031006 | Squared dot spin-up electron-addition matrix element |
| `pair_weight` | 0.0007372261697 | Squared phase-derivative pair matrix element, in squared gap units |
<!-- END generated: static-results -->

At this point the lowest even state is the global ground state. The positive
quantity `signed_gap = E_D - E_S` gives the lowest electron-addition/removal
excitation. When a sweep makes it negative, the doublet becomes the ground state;
the physical excitation energy is then its absolute value. The even-sector pair
excitation is `E_Sstar - E_S` and is evaluated independently of the parity-changing gap.

At zero field the two doublet spin projections are degenerate; the example checks
this. In a field, compute both projections before identifying the global minimum.
For spin-mixing Hamiltonians choose `twice_sz=None` and classify states using the
appropriate observables in the resulting space of all spin projections.

## Interpret observables and array axes

`result.energies[i]`, `result.residuals[i]` and `result.observables[name][i]`
refer to the same zero-based state index. QP vectors are columns of
`result.vectors`, with rows in `result.basis` order. `qp_weights[i, n]` gives the
probability of finding `n` bath quasiparticles in state `i`. Physical dot charge
is a local observable, while total parity labels the full dot-plus-reservoir state
and is specified by the chosen sector.

`bra_result.matrix_elements(operator, ket_result)[i, j]` has bra states as rows
and ket states as columns. Thus the shape of `addition` in the example is one
odd final state by two even initial states. For the pair matrix, the transition
from the lower to the upper even state is `[1, 0]`. Its squared magnitude is a
probe strength. The golden-rule transition rate also depends on the probe's
amplitude and frequency; initial-state populations determine the contributions
to the measured response.

Use named probes such as `d_up_dagger` to obtain physical dot-electron operators
in the Hamiltonian's canonical coordinates. For an arbitrary electron operator,
use `h.physical_operator(...)` with an uncompressed declared electron map.
Multiplying individually projected operators can omit virtual intermediate
states: form the full fermionic product before projecting it.

## Convert the current to physical units

The observable `phase_derivative` and the named probe `current` represent
the derivative of the Hamiltonian along its declared phase coordinate. Multiply
its expectation value by $2e/\hbar$ to obtain the physical current. With all energies
normalized to the gap, the example's number is $I/(2e\Delta/\hbar)$.

Check a branch current against a finite difference of the **same sector's**
energy. Across a ground-state parity crossing, first follow each smooth sector
branch, then take the lower energy. Differentiating an excitation energy gives
a difference of state currents. At the crossing, evaluate the appropriate
one-sided ground-state derivatives, which can differ.

## Sweep and assess the result

![Subgap energies and current of the reference junction](../examples/01_subgap_states/expected.png)

This Tutorial 01 figure uses the same reference parameters and QP cutoff. Its
shading marks the lowest continuum threshold across sectors; determine the threshold
for each probe from its allowed final sector. The full tutorial additionally
checks current derivatives, spin degeneracy, particle–hole symmetry and the
decoupled-dot limit. These are independent physical checks of the selected finite
calculation. Establish suitable cutoffs and bath sizes by refinement in each
parameter regime.

Next use the [scan-and-save workflow](workflows.md), then vary the
[QP cutoff](qp_workflow.md) and [bath resolution](reservoir_workflow.md)
separately. The [result reference](reference/results.md) describes all fields;
the [record reference](reference/records.md) describes their on-disk representation.
