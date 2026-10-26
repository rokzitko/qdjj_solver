# General impurities and custom operators

`reference_model` is convenient for a symmetric one-orbital junction.
`make_model` accepts an explicit impurity, reservoirs and contact matrices,
and supports unequal gaps or couplings, multiple orbitals and spin mixing.
Both return a `Hamiltonian`; choosing a solver remains a separate step.

## Match a Python model to its JSON input

The complete files
[`examples/handbook/custom.py`](../examples/handbook/custom.py) and
[`custom.json`](../examples/handbook/custom.json) describe the same two-orbital
impurity coupled to one superconducting reservoir. Run:

```sh
python -B examples/handbook/custom.py --output runs/custom-python.json
qdjj-solver solve examples/handbook/custom.json --output runs/custom-cli.json --save-states
```

The Python script constructs the operator explicitly, compares it with the JSON
construction, and compares the two sets of calculated energies. The relevant excerpt is:

<!-- BEGIN generated: custom-example (tools/gen_docs.py) -->
<!-- docs-test: skip -->
```python
n = [number(i) for i in range(4)]
impurity = Impurity(2, -.3*(n[0]+n[1]) + .2*(n[2]+n[3]) + n[0]*n[1])
bath = cosh_grid(pairs=1, delta=1., bandwidth=5.)
tunneling = np.array([[.3, 0., .1, 0.], [0., .3, 0., .1]])
h = make_model(impurity, [bath], [tunneling], phases=[0.])
result = solve(h, cutoff=2, sector=Sector(0, 0), options={"eigenpairs": 2})
save_result(result, output, save_states=True)
```
<!-- END generated: custom-example -->

The first orbital has one-electron energy -0.3 and interaction one; the second
has one-electron energy 0.2 and no local interaction. This explicitly written
Hamiltonian retains the specified additive constant. To use the shifted Anderson
energy convention, include its constant explicitly. Supply tunnelling amplitudes in energy units;
convert hybridization widths using each reservoir's normal-state density of states.

## Understand shapes and ordering

For `M` impurity orbitals and `L` reservoirs:

- Impurity modes are `(0 up, 0 down, 1 up, 1 down, ...)`.
- `tunneling[l]` has two reservoir-spin rows and `2*M` impurity-spin columns.
- `phases` has one value in radians per reservoir.
- A direct contact matrix for `(l, m)` has shape `(2, 2)` with `l < m`.
- Reservoir physical electron modes follow the impurity modes, ordered by
  reservoir, signed level and spin.

Unequal lead couplings are specified through the corresponding tunnelling
matrices, using each bath's normal density convention. Off-diagonal spin entries
can break spin-projection conservation, requiring `Sector(parity, None)`.
Request eta or compression reductions only when their stated symmetry conditions
are satisfied. Check those conditions explicitly for a multiorbital device.

## Build operators before projecting them

`create(i)`, `annihilate(i)` and `number(i)` take zero-based mode indices. The
serialized monomials instead use `+(i+1)` for creation and `-(i+1)` for
annihilation. Their products are written in mathematical order, with the
rightmost operator acting first. The algebra normal-orders complete products.

For example, form `number(0)*number(1)` as a fermionic polynomial, then evaluate
its projection. This retains the contributions of virtual intermediate configurations
at the cutoff boundary.

An impurity operator supplied to `Impurity` is in physical impurity coordinates.
`make_model` transforms it along with the Hamiltonian. An arbitrary full-system
electron operator supplied afterwards must be transformed with
`h.physical_operator`. This requires an uncompressed declared electron map.
Built-in observables include the contractions required by any supported
dark-mode restriction.

## Supply general one- and two-body integrals

`Impurity.from_integrals(one_body, interaction, constant)` uses spin-orbital
indices and the convention

```math
H_d=\sum_{ab}h_{ab}d_a^\dagger d_b+
\frac14\sum_{abcd}U_{abcd}d_a^\dagger d_b^\dagger d_d d_c+C.
```

The one-body matrix is square with even dimension; the interaction tensor has
that dimension on all four axes. The factor one quarter assumes antisymmetrized
two-electron integrals. Match operator ordering and normalization before importing
integrals from another program. Small direct-electron-space checks are useful
for validating a new model.

## Direct finite Hamiltonians

`Hamiltonian(operator, nimp, spins, ...)` supports a completely explicit finite
fermionic problem. Declare all modes, identify the first `nimp` as impurity modes,
and retain the desired additive constant. The constructor checks Hermiticity
and parity conservation; the solver checks compatibility with the requested sector.

For later coordinate-sensitive operations, supply consistent coordinate metadata
or retain the model constructors' declarations. Verify that the bases describe
the same physical modes as well as having compatible dimensions.

The [constructor reference](reference/models.md) gives every keyword, and the
[configuration reference](reference/configuration.md) specifies complex arrays,
operator-term records and contact lists. The API reference documents the methods
for conversion, fingerprints and probe construction.
