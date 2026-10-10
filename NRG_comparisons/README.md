# QP, DMRG, and NRG Results

## Model Progression

| Case | Change in physical model | Result and scope |
|---|---|---|
| [Test1](test1/README.md) | One reservoir, zero detuning and field | Doublet ground state; conditional empirical NRG reference; DMRG all-root residual qualification incomplete. |
| [Test2](test2/README.md) | Add detuning $0.5\Delta$, retain zero field | Doublet ground state; both final same-bath DMRG branches qualify; no determined continuum uncertainty. |
| [Test3](test3/README.md) | Add local field $B=0.03\Delta$, retain detuning | Odd-down ground state and four saved subgap states in finite-bath/single-twist pilots; DMRG residuals fail. |
| [Test4](test4/README.md) | Replace one reservoir by two physical L/R leads with $\phi=3\pi/5$; retain detuning and field | Four observed subgap states and branch currents; qualified small-bath QP/DMRG agreement, coarse NRG working means; campaign closed without strict continuum accuracy. |

Test1, Test2, and Test3 use one reservoir; Test4 uses two physical reservoirs.
All cases use $D=100\Delta$, $U=2\Delta$, and total $\Gamma=0.4\Delta$.
For Test4, $\Gamma=\Gamma_L+\Gamma_R$ with $\Gamma_L=\Gamma_R=0.2\Delta$.
The Test4 occupancy parameter $\nu=0.75$ means energy detuning
$U(1-\nu)=0.5\Delta$, not detuning 0.75. Its direct bridges are $W_N=W_S=0$.
Their data and conclusions are independent: **Test1's conditional empirical NRG
uncertainties do not apply to Tests 2, 3, or 4.** Each case's README is its
canonical scientific report; Tests 2-4 put results before `Supplemental material`.

## Energies and Gap

The common signed gap is $g=(E_{\mathrm{odd,lowest}}-E_{\mathrm{even},0})/\Delta$,
equal to $(E_D-E_S)/\Delta$ at zero field. Negative values mean the odd state is
lower. At finite field the odd splitting is reported separately; parity and
total spin projection $M$ do not assert exact SU(2) multiplets.

| Case | Selected QP setting | QP signed $g$ | Selected NRG setting | NRG signed $g$ | QP odd splitting / $\Delta$ |
|---|---|---:|---|---:|---:|
| Test1 | Cosh26, $q=6$ | -0.1743866885 | Extrapolated, 32 twists per $\Lambda$ | -0.1743860 | 0 by symmetry |
| Test2 | Cosh20, $q=7$ | -0.0630007522 | $\Lambda=2$, 8 twists, keep8000 | -0.0629988764 | 0 by symmetry |
| Test3 | Full cosh8 pilot | -0.0759524061 | $\Lambda=2$, $z=1$, keep2000 pilot | -0.0754980794 | 0.0252607480 |
| Test4 | Cosh12 per lead, $q=6$ | -0.1581296660 | $\Lambda=8$, 8 twists, explicitly mixed retention | -0.1580710903 | 0.0248168300 |

These are selected reported approximations, not calculations at equal numerical
accuracy. Cosh sizes count signed normal levels per reservoir; Test4 therefore
has two copies. The upper even excitation is $0.9368718442\Delta$ in Test3's
full-QP pilot and $0.7653557295\Delta$ in Test4's selected QP result. Neither
finite-field study establishes complete continuum subgap coverage.

## Test1 Summary

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

The test1 results concern one reservoir at zero detuning and field. They do not
validate other geometries or establish rigorous continuum error bounds.

## Test2 Summary

The [detuned comparison](test2/README.md) has a doublet ground state with a gap
near $-0.063001\Delta$. On the same cosh-16 bath, DMRG and QP cutoff 8 agree in
branch energies within $2.0\times10^{-14}\Delta$. DMRG residuals pass at singlet
bond dimension 1024 and doublet bond dimension 768. Bath changes still exceed
some requested QP/DMRG observable targets, and the finite-$\Lambda$ NRG means
have no assigned continuum uncertainty.

[Values](test2/output/observables.csv), [comparisons](test2/output/comparisons.csv),
[complete result data](test2/output/results.json), and two figures are provided.
The [presentation generator](test2/README.md#reproduction) requires
no calculation archive or solver execution.

## Test3 Summary

The [finite-field pilot results](test3/README.md) retain both odd projections
$M=\pm1/2$ and **two** even $M=0$ roots, not only the three primary branches.
Same-bath full QP and DMRG agree closely, but all four DMRG states fail residual
qualification at bond 128. Single-$z$ NRG uses a different bath and provides no
continuum reference or error estimate. Only longitudinal spin quantities are
reported. This remains a pilot, not a production convergence study.

## Test4 Summary

The [two-reservoir final report](test4/README.md) uses the canonical physical
left/right model, not stock even/odd channels. The
[model and parameter mapping](test4/README.model_and_parameters.md) specifies
$\phi=\phi_R-\phi_L=3\pi/5$, $\nu=0.75$ (energy detuning $0.5\Delta$),
local field $B=0.03\Delta$, and total $\Gamma=0.4\Delta$.
The baseline has $W_N=W_S=0$; nonzero bridges are validation controls only.
All Test4 NRG calculations, including controls, use `SPU1` and $\Lambda\geq4$.

The final mixed NRG mean uses $\Lambda=8$, `Nmax=16`, `keepenergy=16`, floor400,
and eight equal-weight shared twists: seven keep4000 points and keep5000 only at
$z=0.125$. It is **not a uniform keep5000 grid**. The original uniform keep4000
mean remains separate. Against finite cosh12 QP $q=6$, the four branch currents
in units of $2e\Delta/\hbar$ differ by 0.15-0.42%, both even excitations by about
0.04%, and the odd splitting by 1.08%. These discrepancies are not error bars.

All four cosh4 DMRG states qualify; larger-bath four-state qualification remains
incomplete despite strong same-bath agreement. **The campaign is closed, no
further runs are scheduled, and strict continuum accuracy is not established.**
The older [final-report URL](test4/FINAL_REPORT.md) remains a thin pointer.

## Data and Reproduction

Test1's data and reproduction links above are unchanged. The other cases ship
self-contained scientific datasets and presentation generators:

| Case | Data and tables | Presentation-only reproduction |
|---|---|---|
| Test2 | [JSON](test2/output/results.json), [observables](test2/output/observables.csv), [comparisons](test2/output/comparisons.csv) | [Instructions](test2/README.md#reproduction) |
| Test3 | [JSON](test3/output/results.json), [report](test3/output/report.md), [observables](test3/output/observables.csv), [comparisons](test3/output/comparisons.csv) | [Instructions](test3/README.md#reproduction) |
| Test4 | [JSON](test4/output/results.json), [report](test4/output/report.md), [observables](test4/output/observables.csv), [comparisons](test4/output/comparisons.csv), [twist means](test4/output/z_averages.csv) | [Instructions](test4/README.md#reproduction) |

Presentation regeneration reads shipped data, writes ignored `runs/presentation/`
directories, and needs no solver execution or historical calculation archive.
Each linked reproduction section provides a standard-library-only command.
Fresh numerical calculations are a separate workflow with independent convergence
requirements. Test4's optional manual NRG bundle is exported from canonical
sources, not shipped as a duplicate calculation folder.
