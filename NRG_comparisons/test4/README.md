# Two-Lead Finite-Field Junction: Final Results

## TL;DR

- **Campaign closed:** no further runs. Model/same-bath checks pass, but strict
  continuum accuracy was **not established**, and no continuum error bars are claimed.
- **Four observed subgap states:** odd $M=-1/2$ is lowest, followed by odd
  $M=+1/2$ and two even $M=0$ roots. This is not a completeness certificate.
- **NRG versus QP:** the final mixed-retention NRG mean and finite cosh12 QP
  differ by 0.15-0.42% in currents, about 0.04% in both even excitations,
  and 1.08% in the odd splitting. Neither is an exact continuum reference.
- **QP versus DMRG:** same-bath agreement is strong. All four cosh4 DMRG states
  qualify, but no complete larger-bath four-state DMRG result does.

This follows Test1 (symmetric single lead), Test2 (detuning), and Test3 (local field).
Earlier uncertainties do not transfer; see the [comparison index](../README.md).

## Model

Two identical physical left/right superconducting reservoirs have gap $\Delta=1$,
normal-state **half-bandwidth** $D=100\Delta$, $U=2\Delta$, and **total**
$\Gamma=0.4\Delta$, split as $\Gamma_L=\Gamma_R=0.2\Delta$. The preferred
occupancy $\nu=0.75$ maps to detuning $\epsilon=U(1-\nu)=0.5\Delta$ and
Anderson level $\epsilon_d=-0.5\Delta$. The local field is $B=0.03\Delta$,
$\phi=\phi_R-\phi_L=3\pi/5$, and $W_N=W_S=0$.

```math
H_d=\frac{U}{2}(n_d-1)^2+\epsilon(n_d-1)+B S_d^z,
\qquad S_d^z=(n_\uparrow-n_\downarrow)/2.
```

Positive $B$ raises spin up. Total $M=S_d^z+S_L^z+S_R^z$ is conserved;
finite-field parity/$M$ labels are **not** exact singlet/doublet multiplets.
The [mapping](README.model_and_parameters.md) specifies pairing signs, bridge
controls and the physical L/R [NRG model](model.m).
All NRG calculations use `SPU1` and $\Lambda\geq4$. Cosh bath sizes count signed
normal levels **per lead**; each record is used for both leads. Bath weights
are never renormalized, and `compress=False` retains the full declared bath.

## Energies and Gap

The two NRG means below both use $\Lambda=8$, eight shared lead twists $z=j/8$,
`Nmax=16`, `keepenergy=16`, and `keepmin=400`, with equal weights $1/8$.
**Uniform** means keep4000 at every twist. **Mixed** means seven keep4000 points
and keep5000 **only at $z=0.125$**. It is not a uniform keep5000 grid.
QP uses cosh12 with cutoff $q=6$, another finite approximation.

Energies subtract each method's two disconnected BCS lead vacua and retain the
centered impurity constant. The manuscript convention adds $0.0625\Delta$ to
every energy. Relative absolute-energy discrepancies depend on this reference.

| State | Uniform NRG $E/\Delta$ | Mixed NRG $E/\Delta$ | QP $E/\Delta$ | Mixed/QP discrepancy (%) |
|---|---:|---:|---:|---:|
| Odd down, $M=-1/2$ | -1.098614553 | -1.105084686 | -1.111336132 | 0.566 |
| Odd up, $M=+1/2$ | -1.073621681 | -1.079997766 | -1.086519303 | 0.604 |
| Even root 0, $M=0$ | -0.937495305 | -0.947013595 | -0.953206466 | 0.654 |
| Even root 1, $M=0$ | -0.330967703 | -0.339422916 | -0.345980403 | 1.932 |

| State | Uniform NRG excitation / $\Delta$ | Mixed NRG excitation / $\Delta$ | QP excitation / $\Delta$ | Mixed/QP discrepancy (%) |
|---|---:|---:|---:|---:|
| Odd down | 0 | 0 | 0 | Not defined |
| Odd up | 0.024992871 | 0.025086920 | 0.024816830 | 1.077 |
| Even root 0 | 0.161119248 | 0.158071090 | 0.158129666 | 0.037 |
| Even root 1 | 0.767646850 | 0.765661770 | 0.765355729 | 0.040 |

The signed gap $g=(E_{\mathrm{odd,down}}-E_{\mathrm{even},0})/\Delta$ is
$-0.158071090$ for mixed NRG and $-0.158129666$ for QP. The odd-up excitation
is the field splitting, not this signed gap. Each excitation subtracts the
lowest **saved** energy at its own twist before averaging. Four nondegenerate
subgap states are observed in these inventories, not an imposed three-state count.

Branch currents are $J=I/(2e\Delta/\hbar)=\Delta^{-1}\partial E/\partial\phi$;
using $e\Delta/\hbar$ units would double every number.

| State | Uniform NRG $J$ | Mixed NRG $J$ | QP $J$ | Mixed/QP discrepancy (%) |
|---|---:|---:|---:|---:|
| Odd down | -0.009330124 | -0.009417139 | -0.009397874 | 0.205 |
| Odd up | -0.009627039 | -0.009716108 | -0.009689263 | 0.276 |
| Even root 0 | 0.060131315 | 0.059654293 | 0.059562238 | 0.154 |
| Even root 1 | -0.119825340 | -0.119541100 | -0.119036787 | 0.422 |

Discrepancy (%) is $100\lvert X_{\mathrm{NRG}}-X_{\mathrm{QP}}\rvert/\lvert X_{\mathrm{NRG}}\rvert$,
using **mixed NRG** in the denominator; near zero, use absolute differences.
Digits identify stored values, not achieved precision. Refinement improves the
current and even-excitation comparisons but worsens the odd-splitting discrepancy
from 0.70% to 1.08%. Agreement is not uniformly improved or an error estimate.

## Probabilities and Moment

$P_0$, $P_1$, and $P_2$ are empty, singly occupied (both spins), and doubly
occupied impurity probabilities. $m_d=\langle S_d^z\rangle$ (`m_dot`) is in units of
$\hbar$; the bath moment includes **both complete leads**, $m_b=M-m_d$.

| State | Method | $P_0$ | $P_1$ | $P_2$ | $m_d$ |
|---|---|---:|---:|---:|---:|
| Odd down | Uniform NRG | 0.09614116 | 0.84970237 | 0.05415647 | -0.41433346 |
| Odd down | Mixed NRG | 0.09633208 | 0.84945358 | 0.05421434 | -0.41421988 |
| Odd down | QP | 0.09607791 | 0.84989002 | 0.05403207 | -0.41449547 |
| Odd up | Uniform NRG | 0.09808790 | 0.84712590 | 0.05478620 | 0.41254801 |
| Odd up | Mixed NRG | 0.09832563 | 0.84682361 | 0.05485077 | 0.41240603 |
| Odd up | QP | 0.09801996 | 0.84732164 | 0.05465840 | 0.41271957 |
| Even root 0 | Uniform NRG | 0.59581286 | 0.34769213 | 0.05649501 | -0.00331203 |
| Even root 0 | Mixed NRG | 0.59847227 | 0.34549795 | 0.05602978 | -0.00327380 |
| Even root 0 | QP | 0.59942683 | 0.34460183 | 0.05597134 | -0.00323597 |
| Even root 1 | Uniform NRG | 0.14223785 | 0.57592698 | 0.28183517 | -0.01718569 |
| Even root 1 | Mixed NRG | 0.14198408 | 0.57365714 | 0.28435878 | -0.01702870 |
| Even root 1 | QP | 0.14153456 | 0.57348495 | 0.28498049 | -0.01720049 |

Mixed NRG/QP probability discrepancies are 0.03-0.35%. Odd-state local moments
differ by 0.067-0.076%; the smaller even moments differ by about 1.0-1.2%.
Both methods find predominantly singly occupied odd states, an empty-dominated
lower even state, and a more doubly occupied upper even state. Total $M=0$
does not force the dot moment to vanish: the bath compensates its polarization.
Only longitudinal correlations are reported, formed at each twist before averaging.

## Achieved Accuracy

Targets remain $10^{-9}+10^{-6}\lvert X\rvert$ for QP/DMRG and
$10^{-6}+10^{-3}\lvert X\rvert$ for NRG. Many coarse-NRG averages show subpercent-
to-few-percent **observed stability**, not continuum errors. Twist spread is not an error bar.

Same-bath QP/DMRG agreement validates the finite problem, not bath accuracy.
Larger QP baths retain unresolved cutoff/bath dependence; larger DMRG states
are not all residual-qualified; NRG retains small-twist retention and
finite-$\Lambda$ sensitivity. No strict continuum or spectral-completeness goal
was established. See [Supplemental material](#supplemental-material) for evidence.

## Data

- [Results JSON](output/results.json): normalized scientific dataset (about
  4.65 MiB), including all physical models, bath records, Wilson coefficients,
  numerical specifications, solver checks, outcomes, and uniform/mixed profiles.
- [Generated report](output/report.md): full observables, controls and convergence evidence.
- [Observable table](output/observables.csv), [comparison table](output/comparisons.csv),
  and [twist averages](output/z_averages.csv): stored precision, operands,
  qualifications, discrepancy denominators, and all retained profiles.
- [Physical input](input/physical.json) and [pilot baths](input/baths.json):
  starting definitions. The JSON's `baths` registry contains **all** study baths.

The dataset retains 157 completed records, 2616 original root rows, 2048 additional
ED/QP control roots, and 88 historical noncompletion outcomes (not unique failures).
Public checksums differ from private-original reference hashes. Software versions
and relative source identities remain, not private execution metadata.

## Reproduction

From the repository root, regenerate with Python 3.11+ and only the standard library:

```sh
python -B -S NRG_comparisons/test4/generate.py --tables-only
```

[generate.py](generate.py) reads **only shipped data**, writing JSON, a report and
three CSVs to ignored `runs/presentation/`, or `--output-dir PATH`. It protects
shipped `output/` and rejects output symlinks; `--overwrite` replaces generated
presentations only. It validates public hashes, identities, grids and recorded
qualifications, and replays scalar controls, not raw MPS/HDF5 or solvers.
No historical runs, NRG installation or Mathematica license are needed.

## Supplemental material

### Controls

The 37 controls cover atomic/noninteracting limits, units, spin/phase reversal,
L/R-to-SDS equivalence, both bridge terms, and derivatives of **all** phase links.
Checked finite energies/nondegenerate observables agree among untruncated NRG,
electron ED, full QP and small-system DMRG within $1.1\times10^{-12}$;
phase-derivative discrepancies are below $8.6\times10^{-10}$. At $\Gamma=0$,
a nonzero bridge still carries current; its energy is not subtracted away.
Nonzero bridges are controls only.
The additional $\Lambda=8$ control compares all 1024 energies within
$7.4\times10^{-13}\Delta$; scalar discrepancies are below $4.3\times10^{-12}$.
Degenerate individual expectations remain basis dependent, not precision evidence.

### QP and Bath Refinement

Each four-state refinement checks 14 observables plus excitation, or 60 values.
The following historical finite-solver checks use the QP/DMRG targets above.
The generated comparison CSV preserves those backend-specific requested targets.
Separate 1% and 5% diagnostic tiers use a $10^{-6}$ absolute floor; they do not
replace the tighter finite-solver targets.

| Refinement | Values within finite target / 60 | Largest change / target |
|---|---:|---:|
| Cosh6 QP $q=4$ to 6 | 1 | Above 1 |
| Cosh6 QP $q=6$ to 8 | 60 | 0.151 |
| Cosh6 QP $q=8$ to 10 | 60 | 2.28e-6 |
| Cosh8 QP $q=6$ to 8 | 60 | 0.208 |
| Cosh10 or cosh12 QP $q=4$ to 6 | 1 each | Above 1 |
| Cosh4 full to cosh6 $q=8$ | 1 | 5.14e5 |
| Cosh6 $q=10$ to cosh8 $q=8$ | 1 | 2.17e5 |
| Cosh8 $q=8$ to cosh10 $q=6$; cosh10 to cosh12 $q=6$ | 1 each | Above 1 |

Bath comparisons pass only the identically zero ground excitation. Upper-even
currents are -0.11137839, -0.12830743, -0.12199600, -0.11899745, and -0.11903679
for these cosh4/6/8/10/12 endpoints: convergence is not monotonic. Per-lead weight
sums on cosh2/4/6 are approximately 0.3765/0.9098/0.9949, without normalization.
Small projected residuals do not bound omitted QP sectors; cosh6/8 cutoff
stability cannot be transferred to cosh10/12, and full cosh6 was not calculated.

### DMRG Residuals

Cosh4 at $\chi=256$ gives four qualified states with residuals from
$3.98\times10^{-13}$ to $9.86\times10^{-13}$ and same-bath full-QP agreement
within $8.63\times10^{-7}$ of the finite target. Earlier $\chi=64,128$ states
remain unqualified. The original $\chi=32$ attempt failed to complete and had
canonical-form warnings; its unstarted refinements are not successful results.

On cosh6 at $\chi=1024$, only odd down passes the $10^{-7}\Delta$ threshold:
its residual is $9.851\times10^{-8}\Delta$, versus $1.021668\times10^{-7}\Delta$
for odd up and $1.78150\times10^{-7}\Delta$, $5.43058\times10^{-7}\Delta$ for
the even roots. Rounding is not a pass. Even energies agree with QP within
$2.81\times10^{-13}\Delta$ despite those residual failures. Cosh10 $\chi=512$
and cosh12 $\chi=256$ match QP $q=6$ on all 60 observed targets but remain
residual-unqualified. Resource-limited refinements remain uncompleted.
Qualification also checks norms, orthogonality, sweeps and complete same-call
root prefixes. Separate sectors need compatible Hamiltonians/coordinates and
a qualified ground reference; excited even roots are never spliced across calls.

### NRG Retention and Field Bound

At $W_S=0$, spin reversal gives the necessary exact bound
$\lvert E_{M=+1/2}-E_{M=-1/2}\rvert\leq\lvert B\rvert=0.03\Delta$.
The following remain **truncation-invalid**, even when later refinements pass:

| $\Lambda$ | $z$ | keep | `Nmax` | Odd splitting / $\Delta$ |
|---|---:|---:|---:|---:|
| 4 | 1 | 200 | 20 | 0.052854142 |
| 4 | 1 | 500 | 20 and 22 | 0.037487893 |
| 4 | 0.125 | 3000 | 22 | 0.031688318 |

At $\Lambda=4,z=1$, cap 1500 to 2000 at floor 300 reaches 10.79 current-target
scales; floor 300 to 400 at cap 2000 stays below 0.00459. Shells 2-8 select the
cap, later shells the window/floor, and the last uses `lastall`. Early-shell cap
truncation dominates the tested floor sensitivity. Exact unclamped window counts
and largest individual symmetry blocks are unavailable.

Cap 3000 to 4000 passes all 60 strict targets at $z=1$ (maximum 0.445 target
scales), but only the last two of eight twists pass completely. At $z=0.125,0.25$,
cap 4000 to 5000 passes only 1/60 each. Cap 6000 to 7000 improves to 9/60 and
28/60 strict passes; all pass the diagnostic 5% tier, not an error certificate.
The keep3000 eight-twist mean remains unqualified; averaging cannot repair its bound-failing member.

### Twist, Length and Discretization

At $\Lambda=4$, keep4000 two-to-four and four-to-eight twist means pass only
15/60 and 9/60 strict targets. At $\Lambda=8$, keep3000-to-4000 and four-to-eight
twist means pass 8/60 and 11/60; main non-absolute-energy step sensitivities are
subpercent to about 2.6%, but endpoint current changes reach 19.2% at $z=0.125$.
Retention is better at most twists, not uniformly converged. All tested
$\Lambda=8$ points pass the field bound, which is necessary but not sufficient.

Length probes cover `Nmax=20` to 22 and 22 to 24 at $\Lambda=4$, and 16 to 18
at $\Lambda=8$; window probes cover `keepenergy=14` to 16 and 16 to 18.
Qualified pairs pass all 60 strict checks, including probes at $z=0.125,0.5$.
They do not cure retention error. Cross-$\Lambda$ comparisons change the finite
bath as well as discretization; they establish no continuum extrapolation or
monotonic bias. Lower recorded work at $\Lambda=8$ is not a controlled hardware
benchmark. Incomplete grids are never averaged as if complete.

### Mixed Endpoint and Limitations

The final $\Lambda=8,z=0.125$ cap 4000-to-5000 step passes 1/60 strict, 18/60
1%, and 38/60 5% checks, all with absolute floor $10^{-6}$. Its odd splitting
changes from 0.026244507 to 0.026996893; currents change by up to 7.40%, and
the lower-even excitation by 15.33%. The endpoint is still retention-sensitive.
Replacing its $1/8$ contribution changes mean currents by 0.24-0.93% and even
excitations by 0.26% and 1.93%. All mean-update checks pass 5%, **not a 5%
continuum bound**. Uniform data remain unchanged alongside the mixed profile.

Missing/ambiguous-state flags and saved excluded/unqualified roots remain in the
dataset. The lowest saved energy does not exclude unexamined sectors. Three initially
uncollected NRG results were independently recovered without rerunning; their
original incomplete classifications remain distinct from later verification.
Failed and deferred outcomes do not become convergence evidence. NRG residuals
are unavailable, never zero. The campaign is closed with these limitations.

### Fresh Calculations

Fresh calculations use public baths/specifications and [physics.py](physics.py),
[finite.py](finite.py), [nrg.py](nrg.py), [dmrg_states.py](dmrg_states.py),
[continuum_finite.py](continuum_finite.py) and [retention_diagnostics.py](retention_diagnostics.py).
New states need new qualification.

For a standalone NRG template, export the canonical sources to a **new** folder:

```sh
mkdir -p NRG_comparisons/test4/runs
python -B NRG_comparisons/test4/manual.py export NRG_comparisons/test4/runs/manual-fresh
```

[manual.py](manual.py) exports sources/instructions, not a published run folder.
On Linux, use `looper`, external [NRG Ljubljana](https://github.com/rokzitko/nrgljubljana),
and licensed Mathematica for preparation/execution (not needed for export).
Use `nrginit`, `nrg`, `math` on `PATH`, or `LOOPER`/`NRGINIT`/`NRG`/`MATH` overrides.
Extraction needs NumPy/h5py; QP/DMRG need the solver/optional backend.
Exported `1_run` prepares only; `1_run --execute` also runs; `2_proc` extracts the
complete grid. Saved roots include extras, but means cover only three primary
branches, unlike this four-state report. Defaults (eight twists,
$\Lambda=4$, keep200, floor40, `Nmax=20`, window14) are an **unconverged budget
template**. Review `param.loop` and set limits: there is no time/memory watchdog.
