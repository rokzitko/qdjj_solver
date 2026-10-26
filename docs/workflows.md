# Parameter scans, saved results and reproducible analysis

A reproducible parameter scan retains the results at each point alongside the final curve.
Record the physical model, approximation, sector, diagnostics and failed points
so that a change of ground-state parity can be distinguished from a numerical
failure. This chapter follows the reference junction of
[static calculations](static_calculations.md).

## Run a phase scan

With `.[plots]` installed, run from the source root:

```sh
python -B examples/handbook/scan.py --output runs/handbook/scan
python -B examples/handbook/scan.py --output runs/handbook/scan --plot-only
```

The first command evaluates nine phases between zero and pi, saving both sector
records and a scan table. The second redraws the figure directly from the saved records.
The complete [script](../examples/handbook/scan.py) accepts
`--input` for another reference-junction JSON and `--points` for a denser grid.

The loop below is extracted from that script. `config`, `output` and `points`
are its explicit arguments, and `records` starts as an empty list.

<!-- BEGIN generated: scan-example (tools/gen_docs.py) -->
<!-- docs-test: skip -->
```python
for index, phi in enumerate(np.linspace(0., np.pi, points)):
    point = {**config, "model": {**config["model"], "phi": float(phi)}}
    h = model_from_config(point, default_compress=False)
    row = dict(index=index, phi=float(phi), status="ok")
    try:
        even = solve(h, cutoff=point["cutoff"], sector=Sector(0, 0))
        odd = solve(h, cutoff=point["cutoff"], sector=Sector(1, 1))
        save_result(even, output/f"{index:03d}-even.json")
        save_result(odd, output/f"{index:03d}-odd.json")
        row.update(singlet=float(even.energies[0]), doublet=float(odd.energies[0]),
                   current_singlet=float(np.real(even.observables["phase_derivative"][0])),
                   current_doublet=float(np.real(odd.observables["phase_derivative"][0])),
                   residual=max(float(even.residuals[0]), float(odd.residuals[0])))
    except (ValueError, RuntimeError) as error:
        row.update(status="failed", error=str(error))
    records.append(row)
    write_json(output/"scan.json", dict(input=config, points=records))
```
<!-- END generated: scan-example -->

Results are saved immediately after each completed point. A rejected or failed calculation remains in
`scan.json` together with its error message. Plotting a complete scan requires
successful results at every point. For DMRG, also store and inspect the convergence
flags to identify which states satisfy the finite-problem checks.

## Identify branches before taking a minimum

Follow the lowest even energy and lowest odd energy separately. Their difference
can cross zero; the nonnegative excitation is its absolute value. Construct the
global ground-state current from the lower branch only after retaining both
state-resolved currents. Refine scan density to locate crossings, and refine the
bath and many-body representation to converge the underlying energy branches.

The energy ordering of excited states can change within a sector. To follow a particular state,
compare symmetry labels, observables and overlaps in compatible coordinates.
At exact degeneracy, compare the whole relevant subspace. Use these physical
labels and overlaps to track states as their energy ordering changes.

## Analyse scalar records without a solver

Scalar JSON records suffice for energies, previously calculated observables, diagnostics
and plotting. From the output directory, ordinary Python can read them:

<!-- docs-test: skip -->
```python
import json
from pathlib import Path

record = json.loads(Path("000-even.json").read_text())
print(record["energies"])
print(record["observables"]["phase_derivative"])
print(record["calculation"]["qp_cutoff"])
```

This output-dependent excerpt is tested against a freshly generated scan by the
handbook checks. Resolve paths relative to the chosen output directory.
Complex values use real/imaginary objects;
`qdjj_solver.common.io.complex_array` decodes that representation when needed.

## Reload states when new matrix elements are needed

Wavefunction amplitudes are stored in the companion state file. Reuse the static example's
`even.json`, saved with `save_states=True`, to calculate a new matrix element:

<!-- docs-test: skip -->
```python
from qdjj_solver import load_result

states = load_result("even.json")
charge = states.matrix_elements("charge")
print(abs(charge[1, 0])**2)
```

The companion filename is stored in `state_artifact.file`. Keep it alongside the
JSON when moving the result. QP states use NPZ; DMRG uses TeNPy HDF5. Checksums
associate the companion with the scalar record. Preserve that pair and save
derived results separately to retain a valid checkpoint. `save_result` currently
supports static QP and DMRG eigenstates.

## What belongs in the calculation record?

- The resolved physical model, bath coefficients and energy reference.
- Sector, QP cutoff or MPS settings, eigensolver tolerances and seeds.
- State-resolved energies and observables, with named units in derived tables.
- Residuals and other convergence diagnostics.
- Software/numerical-library information and implementation identity.
- The scan coordinate and a status for every attempted point.

The standard result format already records much of this information. The scan
table adds the relationship between points and preserves failures. A repository
commit identifies the source checkout; the stored implementation digest
identifies the installed solver files actually used. They answer different
provenance questions.

## Continue safely

QP `initial` is a nonzero vector in the target basis ordering. Reusing a vector
after changing cutoff needs an explicit embedding. DMRG `initial` is a prior
result with compatible coordinates and layout. The detailed constraints are in
[starting states and continuation](interfaces.md#starting-states-and-continuation).

A reproducible report also states its limitations: describe the achieved convergence of the
reported quantity and retain unresolved differences. Support continuum-accuracy
estimates with explicit many-body and bath refinements. See the
[configuration](reference/configuration.md) and [record](reference/records.md)
appendices for complete field definitions.
