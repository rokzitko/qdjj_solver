# Finite-Field Superconducting Impurity: Pilot Results

## TL;DR

- **Spin-down ground state:** the full cosh-8 QP pilot has an odd-state splitting
  of $0.0252607480\Delta$ and a lowest even excitation of $0.0759524061\Delta$.
- **Four saved subgap states:** both full cosh-8 QP and single-$z$ NRG find a
  second even state near $0.937\Delta$. It is included below, not discarded to
  enforce the original three-state expectation.
- **QP versus DMRG:** the same-bath energies and observables closely agree, but
  DMRG at bond 128 fails the prescribed residual qualification.
- **Pilot only:** these are finite-bath and single-$z$ outcomes, **not production
  results, a continuum extrapolation, or an uncertainty estimate**. Test1's
  empirical NRG uncertainties do not transfer to this different model.

The progression is Test1 (particle-hole-symmetric single lead), Test2 (detuned,
zero-field single lead), Test3 (detuned, finite-field single lead), then Test4
(two leads with a phase difference).

## Model

One unpolarized superconducting reservoir with gap $\Delta=1$, normal-state
**half-bandwidth** $D=100\Delta$, $U=2\Delta$, total hybridization
$\Gamma=0.4\Delta$, detuning $\epsilon=0.5\Delta$, and local Zeeman energy
$B=0.03\Delta$. The conventional Anderson level is $\epsilon_d=-0.5\Delta$.

```math
H_d=\frac{U}{2}(n_d-1)^2+\epsilon(n_d-1)+B S_d^z,
\qquad \epsilon=\epsilon_d+U/2,
\qquad S_d^z=(n_\uparrow-n_\downarrow)/2.
```

Positive $B$ raises spin up. `even` denotes the lowest even-parity $M=0$ state;
`odd_down` and `odd_up` denote the lowest odd-parity $M=-1/2$ and $M=+1/2$
states, where $M=S_d^z+S_b^z$. These are spin-projection labels, not SU(2)
multiplets. In particular, `even` does **not** mean an exact singlet.

Energies subtract each bath's disconnected BCS vacuum and retain the centered
impurity constant. NRG subtracts its own actual centered-chain vacuum and
detuning to reach the same convention. Excitations subtract the lowest **saved**
energy of each task. See [model conventions](../../docs/models.md).

## Energies and Gap

The best available pilot endpoints are full cosh-8 QP, DMRG on that same bath
at $\chi=128$, and NRG at $\Lambda=2$, $z=1$, `keep=2000`, `Nmax=42`.
Bath sizes count signed normal levels; weights are unrenormalized and
`compress=False`. NRG uses a **different finite bath**, not a twist average.
All table entries are in units of $\Delta$.

| State | Full QP energy | NRG energy | Full QP excitation | DMRG excitation | NRG excitation |
|---|---:|---:|---:|---:|---:|
| `odd_down`, ground | -1.0982848088 | -1.1030100874 | 0 | 0 | 0 |
| `odd_up` | -1.0730240607 | -1.0777685727 | 0.0252607480 | 0.0252607480 | 0.0252415148 |
| Lowest `even` | -1.0223324027 | -1.0275120080 | 0.0759524061 | 0.0759524061 | 0.0754980794 |
| Second even, `sz+0-r1` | -0.1614129646 | -0.1668928310 | 0.9368718442 | 0.9368718442 | 0.9361172564 |

Displayed digits are rounded stored values, not achieved continuum precision.
The signed odd-down minus even gap is the negative of the even excitation;
the odd-state splitting is listed separately because the field lifts the
doublet degeneracy. The fourth root is below the strict subgap threshold
$(E-E_g)/\Delta<1$ in these pilots, but its continuum interpretation is not
established. Close DMRG entries do not override failed residual checks.

## Probabilities and Moment

$P_0$, $P_1$, and $P_2$ are empty, singly occupied (both spins), and doubly
occupied impurity probabilities. $m_d=\langle S_d^z\rangle$ and
$m_b=\langle S_b^z\rangle$ are in units of $\hbar$.

| State | Method | $P_0$ | $P_1$ | $P_2$ | $m_d$ | $m_b$ |
|---|---|---:|---:|---:|---:|---:|
| `odd_down` | Full QP | 0.0906867209 | 0.8580958259 | 0.0512174532 | -0.4217002420 | -0.0782997580 |
| `odd_down` | NRG | 0.0910511326 | 0.8577490712 | 0.0511997962 | -0.4213997381 | -0.0786002619 |
| `odd_up` | Full QP | 0.0923783311 | 0.8558920517 | 0.0517296172 | 0.4203166008 | 0.0796833992 |
| `odd_up` | NRG | 0.0927726927 | 0.8555037052 | 0.0517236022 | 0.4199746398 | 0.0800253602 |
| `even` | Full QP | 0.5399499960 | 0.3692548888 | 0.0907951152 | -0.0033482572 | 0.0033482572 |
| `even` | NRG | 0.5408054631 | 0.3682368687 | 0.0909576682 | -0.0032532738 | 0.0032532738 |
| `sz+0-r1` | Full QP | 0.2395031313 | 0.5258267495 | 0.2346701192 | -0.0316427701 | 0.0316427701 |
| `sz+0-r1` | NRG | 0.2444079866 | 0.5115402783 | 0.2440517351 | -0.0312180043 | 0.0312180043 |

These four states are individually nondegenerate within their saved sectors,
so their scalar diagonal observables are meaningful for the stated finite
problems. The full-QP rows pass finite-problem qualification, including the
second even root. NRG has no residual estimate; its rows are not continuum
reference values. DMRG probabilities and moments are retained in the
[generated report](output/report.md#probabilities-and-moment), explicitly
unqualified rather than used as extra precision evidence.

The lowest even state is locally polarized despite total $M=0$. Same-sector
degeneracies in other saved control/exclusion roots make their individual
expectations basis dependent; the dataset flags them. A tie between different
conserved spin sectors does not create that ambiguity.

## Achieved Accuracy

All 17 control tasks completed and their checks passed. QP cutoff 6 differs
from full cosh-8 QP by at most 0.0285 of the requested target over the three
primary branches and their excitations. This is finite-bath cutoff evidence,
not bath convergence. All four DMRG states displayed above fail residual
qualification at $\chi=128$.

The unchanged targets are $10^{-9}+10^{-6}\lvert X\rvert$ for QP/DMRG and
$10^{-6}+10^{-3}\lvert X\rvert$ for NRG. Small target-normalized differences
are **not error bounds**. No backend has established continuum accuracy,
spectral completeness, or an empirical uncertainty for this finite-field model.

## Data

- [Results JSON](output/results.json): all 27 tasks, 330 saved states, explicit
  request baths and Wilson-chain coefficients, numerical settings, residuals,
  qualifications, 42 control checks, and original source/receipt hashes.
- [Observable table](output/observables.csv): all saved roots and all 12
  longitudinal observables, including flagged and exclusion roots.
- [Comparison table](output/comparisons.csv): 858 primary-observable/excitation
  comparisons retaining both operands, targets, differences and status.
- [Generated report](output/report.md): four-state endpoint comparison first,
  followed by detailed residual, control, refinement and task-status tables.
- [Physical input](input/physical.json), [explicit baths](input/baths.json), and
  [publication provenance](reference/provenance.json): original hashes are
  distinguished from hashes of edited sources and public files.

Private paths, installation locations, execution identities and session records
are not publication data. The original receipts remain separate; their hashes
are not silently relabeled as hashes of the sanitized projection.

## Reproduction

Regenerate all four presentation assets from the **shipped** dataset with only
Python 3.11 or later, from the repository root:

```sh
python -B -S NRG_comparisons/test3/generate.py
```

[generate.py](generate.py) writes to ignored `runs/presentation/`, or to
`--output-dir PATH`. It refuses to overwrite published outputs, including
symlink aliases. No historical run, analyzer dependency chain, solver package,
raw HDF5, NRG or Mathematica installation is required. This regenerates the
presentation, not the numerical calculation.

For **fresh bounded controls and pilots**, the portable [runner](run.py) and
[manual NRG workflow](manual_run/README.md) remain available. On Linux, install
the solver and provide your own NRG Ljubljana and licensed Mathematica tools:

```sh
python -m pip install -e '.[dmrg,dev]'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NRGINIT=/absolute/path/to/nrginit
export NRG=/absolute/path/to/nrg
export MATH=/absolute/path/to/math
python -B NRG_comparisons/test3/run.py controls NRG_comparisons/test3/runs/fresh-pilot
python -B NRG_comparisons/test3/run.py pilot NRG_comparisons/test3/runs/fresh-pilot
python -B NRG_comparisons/test3/analyze.py NRG_comparisons/test3/runs/fresh-pilot \
  --output-dir NRG_comparisons/test3/runs/fresh-report
```

The [raw-receipt analyzer](analyze.py) retains internal integrity checks, but
its writer always passes through the same allowlisted public projection and
privacy guard. Fresh results need their own qualification; the shipped results
do not authorize or imply production calculations.

## Supplemental material

### Refinement Checks

The ten pilot tasks use QP cutoffs 2, 4, 6 and full on `cosh8`; DMRG bonds
32, 64 and 128; and NRG retention 500, 1000 and 2000 at $\Lambda=2$, $z=1$,
`Nmax=42`. NRG also changes `keepmin`; every numerical input is retained.
Neither full cosh8 QP nor highest-retention single-$z$ NRG is a continuum
reference. Larger baths and independent NRG twist/discretization refinements
have not been performed here.

The comparison denominator is the **new** backend's absolute target plus its
relative target times $\lvert X_{\mathrm{new}}\rvert$; ED uses QP targets for
this diagnostic only. Recorded finite-control pass/fail uses a separate
absolute $10^{-8}$ bound. Missing or ambiguous comparisons remain flagged.

### Residuals

At DMRG bond 128, the untruncated finite-Hamiltonian residuals are:

| State | Residual / $\Delta$ | Threshold / $\Delta$ |
|---|---:|---:|
| `even` | 8.9345e-6 | 1e-7 |
| `odd_down` | 1.3994e-6 | 1e-7 |
| `odd_up` | 1.4315e-6 | 1e-7 |
| `sz+0-r1` | 8.2049e-6 | 1e-7 |

Finite-problem qualification is state-specific and includes root/sweep checks.
An unconverged exclusion root does not automatically invalidate a qualified
primary root. QP residuals refer to the projected cutoff Hamiltonian and do
not bound omitted QP sectors. Missing NRG residuals are unavailable, never zero.

### Controls

Actual short Wilson chains are transformed into finite stars for independent
electron-space ED, full QP and DMRG comparisons with SPU1 NRG. Checked energies
and nondegenerate observables agree to better than $2\times10^{-13}$. Atomic
Zeeman sign, spin reversal, zero field, gap/bandwidth units and two-step
Hellmann-Feynman field derivatives pass their separately recorded tolerances.

The runner preserves its bounded sequential supervision: 4 CPUs, 16 GiB
process-tree RSS, 1200 seconds per task, 7200 seconds per campaign, and one
Mathematica instance, plus finite-QP size gates. Controls must pass before
pilots. Changed identities and failed tasks are not silently retried; sibling
campaigns do not reset the authorized budget. Production remains deferred.
Summed saved task time is 1219 seconds and peak sampled RSS is 0.523 GiB;
these pilot measurements are not larger-bath cost estimates.

### Longitudinal Identities

`P_up` and `P_down` are singly occupied spin probabilities, not densities.
The stored observables satisfy

```math
P_0+P_1+P_2=1,\qquad
P_\uparrow=P_1/2+m_d,\qquad P_\downarrow=P_1/2-m_d,
\qquad \langle(S_d^z)^2\rangle=P_1/4,
```

```math
m_b=M-m_d,\qquad
C_{zz}=M m_d-P_1/4,\qquad
\langle(S_b^z)^2\rangle=M^2-2M m_d+P_1/4,
\qquad C_{zz}^{\mathrm{connected}}=m_d^2-P_1/4.
```

Finite backends measure bath operators directly and check these identities;
NRG derives bath moments from fixed-$M$ identities. No transverse correlation,
total-spin squared, or SU(2) factor is inferred.

### Limitations

Extra roots and sectors are exclusion diagnostics, not discarded observations.
The fourth even state needs bath/retention checks before a continuum claim.
Some other NRG roots lie within $10^{-9}\Delta$ of the gap edge, so their
strict-window membership is not robust. The task ground energy is only the
lowest saved energy, not proof that every sector has been excluded.

Full NRG eigenenergy arrays stay in original receipts. Public compact sector
summaries preserve root counts, first omitted/excluded roots and gap-edge
evidence. Raw-artifact hashes are copied, not revalidated during presentation
regeneration. No uncertainty or convergence conclusion is inherited from
Test1, Test2, or another campaign.
