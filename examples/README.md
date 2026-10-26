# Tutorials and examples

This directory holds two kinds of runnable material:

- **Tutorials** (`NN_topic/`): calculations that accompany the worked examples of the
  [physics manual](../docs/manual/README.md). Each has a stand-alone README with the physics in
  brief, the commands, the expected results and suggestions for variations. Each has a
  `run.py` with a quick profile (seconds) and a full profile (the complete calculation).
  The expected results are checked by `python tools/check_tutorials.py`.
- **Interface examples** (`*.py`, `*.json`): small scripts and input files that show how to
  set up and check calculations with the Python interface and the command line.

Install from source following the root README's
[installation prerequisites](../README.md#install),
including a C++ compiler and the Python development headers. Run all
commands from the main source directory.

## Tutorials

<!-- BEGIN generated: tutorials (tools/gen_docs.py) -->
| Tutorial | What it shows | Physics manual | Quick profile |
|---|---|---|---|
| [00: A first calculation](00_first_calculation/README.md) | Lowest odd-parity state of a small junction with the QP expansion, DMRG and exact diagonalization. | `sec:first-calculation`, `sec:qp` | about 1 s |
| [01: Subgap states, singlet and doublet](01_subgap_states/README.md) | Singlet, excited singlet and doublet of the reference junction versus phase and detuning, the singlet–doublet transition and the supercurrents. | `sec:subgap`, `ch:model`, `ch:static` | about 4 s |
| [02: The gate-driven isolated dot](02_driven_dot/README.md) | Quasienergies, Bessel-function harmonic weights in two gauges, replicas and Floquet zones, and Tien–Gordon sidebands of a gate-driven dot, checked against the exact solution. | `ch:ex-dot`, `sec:gauge`, `ch:probes` | about 2 s |
| [03: The proximitized dot as a driven two-level system](03_two_level_system/README.md) | The Autler–Townes avoided crossing compared with the rotating-wave model, the Bloch–Siegert shift and Bessel-modulated multiphoton resonances of the gate-driven proximitized dot. | `ch:ex-tls`, `sec:driven-largegap`, `sec:rwa` | about 6 s |
| [04: The high-frequency limit](04_high_frequency/README.md) | Coherent destruction of tunnelling, a drive-induced singlet–doublet crossing and the accuracy of van Vleck effective Hamiltonians for a fast phase drive. | `ch:ex-highfreq`, `sec:van-vleck`, `sec:jacobi-anger` | about 13 s |
| [05: Numerical convergence of a driven calculation](05_convergence/README.md) | The Floquet cutoff and its residual bound, a cross-check with the one-period propagator, the static basis against the full Fock space, and bath-size effects at the one-photon threshold of the continuum. | `ch:ex-numerics`, `sec:sambe-method`, `sec:propagator`, `sec:continuum` | about 7 s |
| [06: The phase-driven junction](06_phase_drive/README.md) | Rabi splitting of the pair transition under a phase drive compared with the rotating-wave model, the AC Stark shift of the resonance, the eta selection rule and the Floquet spectrum versus drive amplitude. | `ch:ex-phase`, `sec:jacobi-anger`, `sec:selection`, `sec:rwa`, `sec:regimes` | about 14 s |
| [07: The gate-driven junction](07_gate_drive/README.md) | Rabi splitting of the gate-driven pair transition compared with a phase drive, the complementary eta selection rules versus the operating point, the AC Stark shift and one- and two-photon resonances. | `ch:ex-gate`, `sec:selection`, `sec:regimes`, `sec:rwa` | about 15 s |
| [08: Two-tone spectroscopy of the pair transition](08_two_tone/README.md) | Probe lines of the dressed pair transition versus pump detuning, the Mollow pattern and the Autler–Townes doublet of a gate-pumped junction. | `ch:ex-twotone`, `ch:probes`, `sec:mollow` | about 3 s |
| [09: Photon-assisted tunnelling into a driven junction](09_photon_assisted_tunnelling/README.md) | Photon sidebands of the tunnelling line of a gate-driven junction compared with the simple Tien–Gordon formula and its adiabatic generalization, including the asymmetry of the sidebands $q=\pm1$. | `ch:ex-pat`, `ch:probes`, `ch:ex-dot` | about 3 s |
| [10: Electron spin resonance of the doublet](10_spin_resonance/README.md) | Effective g-factor, Rabi splitting and Bloch–Siegert shift of the Zeeman-split doublet under a transverse drive, and the Mollow and Rabi lines of weak spin probes. | `ch:ex-esr`, `sec:esr-probe`, `sec:mollow` | about 8 s |
| [11: DMRG bond-dimension convergence](11_dmrg_convergence/README.md) | Refinement of an initially unconverged MPS by increasing the bond dimension, with convergence diagnostics, ED comparisons for the same bath and continuation from saved states. | `ch:practical`, `sec:dmrg` | about 8 s |
<!-- END generated: tutorials -->

## Interface examples

### Base-only examples

These commands need only the base installation (`python -m pip install .`);
TeNPy and HDF5 are not required:

```sh
python -B examples/quickstart.py
python -B examples/sector_convergence.py
python -B examples/qp_convergence.py
python -B examples/floquet.py
qdjj-solver solve examples/qp.json --output results/qp.json --save-states
```

`quickstart.py` checks the README's finite-model energy, residual, and
normalization. `sector_convergence.py` compares both parities, refines the QP
cutoff, and computes an excitation energy and current; see its
[worked interpretation](../docs/numerics.md#worked-parity-and-current-example).
`qp_convergence.py` follows a spin-up doublet through increasing QP cutoffs up to
full finite-bath ED, reporting energy, local spin, dot charge probabilities, and a
direct dot–bath spin correlation. It checks the SU(2) doublet identity at every
cutoff. The [cutoff convergence guide](../docs/qp_convergence.md) adds a larger
reference study.
`floquet.py` drives the junction phase at the lowest even-parity resonance and
compares Sambe diagonalization in the Fock and static bases, the one-period
propagator and the two-level rotating-wave model. Its self-checks compare the
splitting with the weak-drive estimate $A\lvert\langle 0\vert\partial_\phi H\vert 1\rangle\rvert$.
See [Floquet calculations](../docs/floquet.md).

### DMRG-required examples

Before running these commands, install the `[dmrg]` extra with
`python -m pip install '.[dmrg]'`. It supplies TeNPy and HDF5 support and is
required for both the DMRG commands and `multi_orbital.py`:

```sh
qdjj-solver solve examples/dmrg.json --output results/dmrg.json --save-states
qdjj-solver solve examples/dmrg_refined.json --resume results/dmrg.json --output results/refined.json --save-states
python -B examples/multi_orbital.py
```

`qp.json` and `dmrg.json` use identical physical coordinates. The QP
example imposes a finite cutoff; change it to `null` for unrestricted finite
space. The multiorbital example compares unrestricted ED and DMRG for a small
spin-mixing model. These examples show how to set up and check a calculation;
their small baths serve as finite-model illustrations.

`dmrg_refined.json` uses the same physical model and MPS site ordering as
`dmrg.json`, raises the maximum bond dimension to 64, and tightens the numerical
tolerances. The `--resume` command
initializes a new optimization from the first run's saved MPS states. Both
commands save wavefunctions so either completed result can be loaded again.

These files are supplied with the source code. To run outside the source tree,
the [main README](../README.md#command-line)
provides a complete JSON input that can be saved in any working directory.
