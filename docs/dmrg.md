# DMRG states and bond-dimension refinement

DMRG represents a correlated state as a matrix-product state (MPS). The bond
dimension limits the entanglement represented by the MPS while retaining all bath-QP numbers.
The same finite `Hamiltonian` can be passed to QP and DMRG. Matching its bath,
coordinates, compression and sector is essential for a meaningful comparison.

## Start with a small reference problem

Install `.[dmrg,plots]` from the source root and run:

```sh
python -B examples/11_dmrg_convergence/run.py --output runs/dmrg-ladder
```

This [complete tutorial](../examples/11_dmrg_convergence/README.md) uses the
small installation-test Hamiltonian, retains all bath modes, and compares with
full ED of the same finite Hamiltonian. It uses a sequence of increasing bond limits, saves each completed
step, and finally saves, reloads and continues the MPS. The full profile adds
intermediate and final steps:

```sh
python -B examples/11_dmrg_convergence/run.py --profile full --output runs/dmrg-ladder-full
```

Both commands use a fixed seed and one computational thread. Small bond limits
are intentionally inadequate: TeNPy can print truncation/canonical-form notices
at those steps. The saved diagnostics show how the approximation improves. The final
result must pass the residual and empirical-stability checks at the stated
tolerances. All diagnostic messages remain visible throughout the refinement.

## Interpret the refinement sequence

![Energy and residual versus MPS bond limit](../examples/11_dmrg_convergence/expected.png)

The energy discrepancy is relative to full ED on the same bath. The residual
is the norm of the full finite-Hamiltonian eigenvalue equation. Values below
1e-14 are placed at that plotting floor. The two curves quantify different aspects
of finite-Hamiltonian accuracy; bath refinement assesses continuum accuracy.

The first coarse state is inaccurate and marked unconverged; the refined energy
agrees with ED within the checked 1e-9 tolerance. `history.json` records the
requested and actual bond dimensions, energies, charge, residuals and overlap errors for
each completed step. The final assessment is also saved in `results.json`.

## Which controls should be changed?

| Question | Controls and evidence |
|---|---|
| Is the state representation large enough? | Increase `chi_max`; compare energy and the observable being reported |
| Has each optimization converged? | Inspect sweep histories; vary the number of sweeps and local-solver precision |
| Is an excited state reliable? | Use independent trials; inspect energy ordering and the full Gram matrix |
| Is the finite eigenvalue equation satisfied? | Evaluate `residuals`; inspect `finite_problem_converged` |
| Does the bath represent the target continuum? | Refine the bath separately, with fixed physical parameters |

`require_convergence=False`, the default, returns coarse states with their
diagnostics. Use `True` for a calculation that must satisfy finite-problem
criteria. The residual calculation has temporary memory costs and can be more
demanding than an expectation value. Turning it off leaves convergence unverified.

The local Lanczos probability tolerance and the full-state residual tolerance
are distinct. Convergence requires both adequate bond dimension and sufficiently
accurate local solves. Increasing the bond dimension can reveal an optimization precision
limit that was previously hidden by truncation.

## Use the existing refinement helper

The tested tutorial uses the following code. Here `h` and `sector` specify the
finite problem, and the function `archive`, defined in the complete script, saves each completed step:

<!-- BEGIN generated: dmrg-example (tools/gen_docs.py) -->
<!-- docs-test: skip -->
```python
options = SolverOptions(eigenpairs=1, seed_trials=1, seed=1729,
                        max_sweeps=40, min_sweeps=8, mixer_sweeps=4,
                        lanczos_probability_tolerance=1e-22)
chis = [2, 4, 16, 32] if profile == "quick" else [2, 4, 8, 16, 32, 64]
result, assessment = converge_bond_dimension(
    h, chis, sector=sector, options=options, stable_steps=1 if profile == "quick" else 2,
    energy_tolerance=1e-9, observable_tolerance=1e-8, callback=archive)
if not assessment["empirical_convergence"]:
    raise RuntimeError("bond ladder did not converge; inspect history.json")
save_result(result, output/"state.json", save_states=True)
restored = load_result(output/"state.json")
continued = solve(h, sector=sector, options=SolverOptions(
    eigenpairs=1, chi_max=chis[-1], seed_trials=1, require_convergence=True,
    lanczos_probability_tolerance=1e-22), initial=restored)
```
<!-- END generated: dmrg-example -->

`converge_bond_dimension` runs the supplied sequence of increasing bond limits and returns its
**last result**, together with an assessment. It deliberately permits coarse
unconverged steps, then checks finite-problem convergence and successive changes
in energies and measured observables. The `empirical_convergence` flag reports
whether the sequence met these criteria. Two stable steps are normally required;
the quick teaching profile requests one, while the full profile requests two.

## Request excited states deliberately

Set `eigenpairs` to the number of low states needed in the selected sector.
They are optimized sequentially with penalties against previously found states.
Different initial trials can converge to different candidate eigenstates. Inspect the
diagnostics for each state and the overlap matrix, especially near degeneracies. Comparing
retained subspaces is more meaningful than comparing individual eigenvectors
inside an exactly degenerate multiplet.

`compare_states(previous, current)` reports a maximum-overlap assignment,
energy/observable changes and singular values of the overlap matrix. Those
singular values are invariant under rotations within the compared retained
subspaces. Establishing the completeness of the low-energy spectrum requires
additional state searches or comparison with an independent reference.

## Save and continue

`save_result(..., save_states=True)` creates JSON plus a named HDF5 companion.
`load_result` checks the record and restores the states; pass that result as
`initial` to a new DMRG solve. The new call supplies its own tolerances and bond
limit. Compatible continuation requires matching canonical coordinates, mode
ordering, grouping, sector and types of conserved quantum numbers.

A checkpoint supplies completed states as initial wavefunctions for a new
optimization. Bath changes usually change the canonical coordinates, so perform a
new optimization when refining the bath. A parameter change is allowed when
those coordinates and the MPS layout remain identical.

The [numerical guide](numerics.md#residuals-and-convergence-checks) explains the
diagnostic definitions, the [option reference](reference/options.md#dmrg-solver)
lists every control, and the [record reference](reference/records.md) specifies
saved fields and companion artifacts.
