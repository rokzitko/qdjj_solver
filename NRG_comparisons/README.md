# QP, DMRG, and NRG Results

## TL;DR

For a superconducting Anderson impurity, the two larger-bath QP calculations
agree closely with each other and with NRG. DMRG reproduces QP ground-branch
results on the identical smaller bath, but that bath differs noticeably from
the larger-bath results.

**NRG meets its conditional empirical accuracy targets. The stricter QP/DMRG
continuum accuracy target remains unestablished.**

Read the [final results, tables, and plots](test1/README.md).
Machine-readable [values](test1/output/observables.csv),
[comparisons](test1/output/comparisons.csv), and
[complete result data](test1/output/results.json) are included.
The [reproduction instructions](test1/README.md#reproduction) retain the parameter
inputs, solver/analysis scripts, and table/figure generator. Fresh calculations
write to ignored `runs/` directories; raw outputs and checkpoints are not shipped.

These results concern one reservoir at zero detuning and field. They do not
validate other geometries or establish rigorous continuum error bounds.
