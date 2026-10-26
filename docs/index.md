# qdjj_solver documentation

`qdjj_solver` computes the subgap states of superconducting quantum-impurity models, such as
quantum-dot Josephson junctions, with a truncated quasiparticle expansion or with DMRG, for
static and periodically driven Hamiltonians.

The documentation consists of three parts, each addressing a different question. Where topics overlap,
one document gives the authoritative account and the others link to it.

| | Answers | Contents |
|---|---|---|
| **This manual** (user and developer manual) | *How do I run, check, save and extend the code?* | Installation, workflows, Python API, command line and input files, results and persistence, convergence practice, troubleshooting, reference tables, developer guide |
| **[Physics manual](manual/README.md)** (LaTeX, PDF) | *What is computed and what does it mean?* | Models and equations, methods and their error control, physical interpretation, worked physics examples with figures |
| **[Tutorials](../examples/README.md)** (`examples/`) | *Show me, and let me run it.* | Runnable, checked calculations with expected results; each README points back to the physics manual |

The [research applications](../applications/README.md) reproduce published results and
document their parameter mappings, references and convergence evidence.

**PDF editions:** releases built with the two-manual workflow include the separate
User's Guide and Reference Manual and physics manual. See [reading and building
the PDFs](user_manual.md), or start with the [first checked calculation](getting_started.md).

```{toctree}
:caption: Getting started
:maxdepth: 1

Installation and first calculation </README>
/docs/getting_started
```

```{toctree}
:caption: User guide
:maxdepth: 2

/docs/models
/docs/interfaces
/docs/static_calculations
/docs/reservoir_workflow
/docs/qp_workflow
/docs/dmrg
/docs/workflows
/docs/custom_models
/docs/floquet
/docs/spectroscopy
/docs/bath_representations
/docs/chain_expansion
/docs/numerics
/docs/qp_convergence
/docs/troubleshooting
```

<!-- BEGIN generated: tutorial-toctree (tools/gen_docs.py) -->
```{toctree}
:caption: Tutorials
:maxdepth: 1

Overview </examples/README>
00: A first calculation </examples/00_first_calculation/README>
01: Subgap states, singlet and doublet </examples/01_subgap_states/README>
02: The gate-driven isolated dot </examples/02_driven_dot/README>
03: The proximitized dot as a driven two-level system </examples/03_two_level_system/README>
04: The high-frequency limit </examples/04_high_frequency/README>
05: Numerical convergence of a driven calculation </examples/05_convergence/README>
06: The phase-driven junction </examples/06_phase_drive/README>
07: The gate-driven junction </examples/07_gate_drive/README>
08: Two-tone spectroscopy of the pair transition </examples/08_two_tone/README>
09: Photon-assisted tunnelling into a driven junction </examples/09_photon_assisted_tunnelling/README>
10: Electron spin resonance of the doublet </examples/10_spin_resonance/README>
11: DMRG bond-dimension convergence </examples/11_dmrg_convergence/README>
```
<!-- END generated: tutorial-toctree -->

```{toctree}
:caption: Research applications
:maxdepth: 1

/applications/README
/applications/bargerbos_2022_parity_diagram/README
/applications/bobok_2025_chain_expansion/README
/applications/choi_2004_kondo_josephson/README
/applications/hecht_2008_gap_edge/README
/applications/kadlecova_2017_asymmetry/README
/applications/paaske_2023_surrogate_spectrum/README
/applications/zalom_2024_multiterminal/README
/applications/zitko_2023_knight_shift/README
/applications/zonda_2016_unequal_gaps/README
/applications/zonda_2023_double_dot/README
/large_Gamma/README
/NRG_comparisons/README
```

```{toctree}
:caption: Reference
:maxdepth: 1

/docs/reference/api
/docs/reference/options
/docs/reference/models
/docs/reference/results
/docs/reference/cli
/docs/reference/observables
/docs/reference/conventions
/docs/reference/concordance
/docs/reference/configuration
/docs/reference/records
```

```{toctree}
:caption: Developer guide
:maxdepth: 1

/CONTRIBUTING
/docs/developer
/docs/validation
/docs/documentation
/docs/releases
```

```{toctree}
:caption: Physics manual and citation
:maxdepth: 1

/docs/manual/README
/docs/user_manual
/docs/references
```

```{toctree}
:hidden:

/docs/benchmarks/baths/README
/docs/benchmarks/baths/tables
/docs/benchmarks/qp_convergence/README
/docs/benchmarks/qp_convergence/table
/NRG_comparisons/test1/README
/NRG_comparisons/test2/README
/NRG_comparisons/test2/manual_run/README
/NRG_comparisons/test3/README
/NRG_comparisons/test3/manual_run/README
/NRG_comparisons/test3/output/report
/NRG_comparisons/test4/README
/NRG_comparisons/test4/README.model_and_parameters
/NRG_comparisons/test4/FINAL_REPORT
/NRG_comparisons/test4/output/report
/large_Gamma/output/report
```
