# Manual Test2 NRG Sweep

This folder is standalone: copy it to a calculation directory if desired. It
does not import `qdjj_solver`, use the test1 runner, or read files outside this
folder.

## Run

Requirements: Linux, Python 3.11 or newer, Perl and `looper`, NRG Ljubljana with
SPSU2/HDF5 support, and Mathematica. Processing additionally needs NumPy and
h5py (`python3 -m pip install numpy h5py`). The report used NRG Ljubljana
2026.09 and Mathematica 13.3.
Load your usual calculation environment before running these scripts.

From this folder:

```sh
./1_run
./2_proc
```

`1_run` first uses `looper` to prepare `1/param`, ..., `8/param`, validates the
complete grid, and then runs `nrginit` and `nrg` sequentially in those directories.
It stops at the first failed calculation. `2_proc` reads every completed result,
prints the averages alongside the report numbers and their signed differences,
and writes the tables described below. It never launches a solver.

For separate preparation and manually selected calculations:

```sh
./1_run --prepare-only
```

Then, from inside each numbered directory, run:

```sh
../run_one
```

Run `./2_proc` from the parent folder after all twists have completed. The
wrappers `1_run` and `2_proc` also work when invoked by path from another working
directory; `run_one` deliberately operates in the current numbered directory.
Do not execute `param.loop` directly: `1_run` checks numbering and saves the
input manifest that `run_one` and `2_proc` require.

The tools are resolved from `PATH`. Optional environment variables `LOOPER`,
`NRGINIT`, `NRG`, `MATH`, and `PYTHON` select individual executable paths, not
shell command strings. No site-specific module is loaded automatically. For example:

```sh
NRG_THREADS=4 ./1_run
```

The default is four MKL/OpenMP threads per NRG calculation, with OpenBLAS 
limited to one. The initializer receives an explicit
Mathematica executable and a per-directory `scratch/` via `NRG_WORKDIR` and
`TMPDIR`. Do not launch multiple initializers unless your license and machine
resources permit it. Logs are `nrginit.log` and `nrg.log` in each run directory.

## Inputs

Edit `param.loop` **before** preparing a fresh sweep. Change `$Nz = 8` in its
prelude to select the number of twists. `Nz` inside each generated `param` is
derived from this one assignment. The integer loop gives

```math
z_i=i/N_z,\qquad i=1,\ldots,N_z.
```

Thus `1/` is $z=0.125$, `2/` is $z=0.25$, and `8/` is $z=1$ by default. This is
not a midpoint grid. The only swept parameter is $z$; all directories share one
value of `Lambda`.

| Setting | Default |
|---|---:|
| `Lambda` | 2 |
| `Nz` | 8 |
| `keep` | 8000 |
| `keepenergy` | 14 |
| `keepmin` | 1000 |
| `Nmax` | 42 |

`Lambda`, `Nz`, `Nmax`, `keep`, `keepmin`, and `keepenergy` are
editable. Choose `Nmax` explicitly when changing `Lambda`; no chain-convergence
claim is inferred. The report used `Nmax=48` at `Lambda=1.8` and `Nmax=34` at
`Lambda=2.5`. This folder bundles reference values only for `Lambda=2`,
`keep=8000`, and four/eight twists. Other numerical settings still produce
averages, but unmatched settings have blank report-comparison columns.

The extractor deliberately fixes test2 physics: one flat-band superconducting
reservoir, $\Delta=1$, half-bandwidth $D=100$, $U=2$, total $\Gamma=0.4$,
detuning $\epsilon=0.5$, and zero field. In the input, `bcsgap=0.01` is
$\Delta/D$ and `delta=0.5` is detuning, **not** the conventional Anderson level
$\epsilon_d=-0.5$. Extending this folder to another physical model requires
updating the extractor and reference data, not just relabeling the results.

## Results and Conventions

`processed/per_z.csv` contains all ten observables for each numbered run.
`processed/averages.csv` contains the means, report values, signed differences
`fresh - report`, and absolute differences divided by the requested NRG tolerance
$10^{-6}+10^{-3}\lvert X_{\mathrm{fresh}}\rvert$.
`processed/summary.json` retains full-precision values, settings, selected
directories, extraction diagnostics, artifact hashes, tool paths, timestamps,
and the reference source. There are no estimates of continuum uncertainty.

Every mean gives equal weight to the complete ordered grid. With eight twists,
processing also averages directories `2,4,6,8` to reproduce the report's four-twist
point. More generally, a sweep whose size is a larger multiple of four also
provides this nested four-twist mean.

The ten observables are `signed_gap`, `singlet.energy`, `doublet.energy`,
`singlet.P0/P1/P2`, `doublet.P0/P1/P2`, and `doublet.moment`. They are extracted
from the lowest singlet and doublet **eigenstates** in the final `raw.h5` shell,
not from the thermal expectations in `custom` or `customfdm`.

`P1` includes both spin orientations. For each branch,

```math
P_0=1-\langle n_d\rangle+\langle n_{d\uparrow}n_{d\downarrow}\rangle,
\qquad P_1=\langle n_d\rangle-2\langle n_{d\uparrow}n_{d\downarrow}\rangle,
\qquad P_2=\langle n_{d\uparrow}n_{d\downarrow}\rangle.
```

`doublet.moment` is $\langle S_d^z\rangle$ in the total spin-up doublet, computed
as `sigma_d_reduced/sqrt(3)`. An isolated up electron gives $1/2$. The singlet
moment is zero by symmetry and is not a separate accuracy target. The template
retains `spect=sigma_d-sigma_d` because listing `sigma_d` in `ops` alone does not
enable the required triplet-operator recalculation.

Absolute energies require a bath subtraction for **each twist separately**.
The normal Wilson chain is reconstructed from `data`, excluding the final unused
hopping coefficient. For its normal eigenvalues $\xi_a$, all in gap units,

```math
E_{\mathrm{vac}}(z)=-\sum_a\sqrt{\xi_a(z)^2+1},
\qquad E_b(z)=E_b^{\mathrm{NRG,total}}(z)-E_{\mathrm{vac}}(z)-0.5.
```

The centered NRG bath convention needs no additional sum of normal eigenvalues.
The last $-0.5$ maps the impurity Hamiltonian to the report's energy reference.
No separate CLEAN calculation is necessary for this analytic subtraction.
`signed_gap` is $(E_D-E_S)/\Delta$; negative means the doublet is lower.

The main eight-twist report values are:

```text
singlet.energy  = -1.022548453165136
doublet.energy  = -1.0855473295931404
signed_gap     = -0.06299887642800428
```

The full reference snapshot is in `reference.json`. Stored digits and small
comparison ratios are not error bounds. Flow completion, retention convergence,
twist convergence, and the continuum limit are distinct checks.
