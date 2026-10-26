---
orphan: true
---

# User's Guide and Reference Manual

This handbook explains how to set up calculations, solve models, assess numerical accuracy and save results with
`qdjj_solver`. It assumes basic Python and condensed-matter physics. The companion
physics manual develops the Hamiltonians, subgap-state physics and driven systems.

Read the user guide in sequence for a first calculation, then use the reference
chapters to look up parameters, return values and input formats. To obtain continuum
results at the desired accuracy, check convergence of the eigensolver, the many-body
representation and the bath discretization separately.

The HTML edition and this PDF use the same canonical chapters. Source-file links
identify the corresponding software revision; the matching source distribution
contains the examples and public inputs needed to run them offline.

```{toctree}
:maxdepth: 1
:caption: Contents

/docs/getting_started
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
/docs/troubleshooting
```

```{toctree}
:maxdepth: 1
:caption: Reference

/docs/reference/options
/docs/reference/models
/docs/reference/results
/docs/reference/observables
/docs/reference/cli
/docs/reference/conventions
/docs/reference/concordance
/docs/reference/configuration
/docs/reference/records
/docs/reference/api
/docs/references
```
