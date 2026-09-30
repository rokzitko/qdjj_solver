# Citation and method references

## Software citation

Please identify the software version used. A ready-to-copy citation for the
current version is:

> Teodor Iličin and Rok Žitko. *qdjj_solver: QP and DMRG solvers for superconducting
> quantum impurities*, version 0.2.0. https://github.com/rokzitko/qdjj_solver.

For BibLaTeX:

```bibtex
@software{qdjj_solver_020,
  author  = {Iličin, Teodor and Žitko, Rok},
  title   = {{qdjj\_solver}: {QP} and {DMRG} solvers for superconducting quantum impurities},
  version = {0.2.0},
  url     = {https://github.com/rokzitko/qdjj_solver}
}
```

[CITATION.cff](../CITATION.cff) contains the citation information in a format
readable by reference-management tools. For modified or unreleased source code,
also record its Git revision identifier (`git rev-parse HEAD`). Saved results
include file checksums identifying the solver code actually used. No
archival DOI is currently assigned to this software.

## Which methods to cite

- **QP expansion:** cite the [QP variational method](#qp-variational-method) as
  well as the software version. State the QP cutoff, symmetry sector, finite bath,
  and canonical bath reference;
  [numerics](numerics.md#qp-diagonalization) defines the projection and the
  required convergence checks.
- **Fitted surrogate reservoirs:** cite Baran, Frost, and Paaske when using
  `fit_surrogate`. Report the fitting mesh and criterion alongside the bath size.
- **Padé chain-expansion reservoirs:** cite Bobok, Frk, Pokorný and Žonda when
  using `chain_expansion`. Report the chain length and finite-/wide-band target.
- **DMRG:** cite the DMRG method references and TeNPy below when using this solver.
  Record the bond dimension, eigenstate-search settings, and convergence evidence.

The references below describe those numerical components; the software version
specifies the implementation used in a calculation.

## QP variational method

- Teodor Iličin and Rok Žitko, *Quasiparticle-resolved variational theory of
  Andreev spin qubits*, manuscript. This is the ASQ method manuscript associated
  with the [archived reference-junction calculations](qp_convergence.md).
  No public DOI or arXiv identifier is recorded here.
- Teodor Iličin and Rok Žitko, *Effects of electron-electron interaction and
  spin-orbit coupling on Andreev pair qubits in quantum dot Josephson junctions*,
  SciPost Physics **20**, 147 (2026).
  DOI: [10.21468/SciPostPhys.20.5.147](https://doi.org/10.21468/SciPostPhys.20.5.147);
  [arXiv:2512.23015](https://arxiv.org/abs/2512.23015).
  Section 3.3 and Appendix B extend the two-QP variational approach to a junction.
- Teodor Iličin and Rok Žitko, *Variational solution of the superconducting
  Anderson impurity model and the band-edge singularity phenomena*,
  SciPost Physics **19**, 006 (2025).
  DOI: [10.21468/SciPostPhys.19.1.006](https://doi.org/10.21468/SciPostPhys.19.1.006);
  [arXiv:2503.18902](https://arxiv.org/abs/2503.18902).
  This is the earlier single-reservoir two-QP construction.

## Method references

1. V. V. Baran, E. J. P. Frost, and J. Paaske, *Surrogate model solver for
   impurity-induced superconducting subgap states*, Physical Review B **108**,
   L220506 (2023). DOI: [10.1103/PhysRevB.108.L220506](https://doi.org/10.1103/PhysRevB.108.L220506).
2. S. R. White, *Density matrix formulation for quantum renormalization groups*,
   Physical Review Letters **69**, 2863 (1992).
   DOI: [10.1103/PhysRevLett.69.2863](https://doi.org/10.1103/PhysRevLett.69.2863).
3. U. Schollwöck, *The density-matrix renormalization group in the age of matrix
   product states*, Annals of Physics **326**, 96–192 (2011).
   DOI: [10.1016/j.aop.2010.09.012](https://doi.org/10.1016/j.aop.2010.09.012).
4. J. Hauschild and F. Pollmann, *Efficient numerical simulations with Tensor
   Networks: Tensor Network Python (TeNPy)*, SciPost Physics Lecture Notes **5**
   (2018). DOI: [10.21468/SciPostPhysLectNotes.5](https://doi.org/10.21468/SciPostPhysLectNotes.5).
5. D. Bobok, L. Frk, V. Pokorný and M. Žonda, *Scalable effective models for
   superconducting nanostructures: Applications to double, triple, and quadruple
   quantum dots*, Physical Review B **112**, 205418 (2025).
   DOI: [10.1103/mxsl-fc96](https://doi.org/10.1103/mxsl-fc96).

TeNPy supplies MPS contractions and the optional DMRG engine. Consult its
[documentation](https://tenpy.readthedocs.io/) for current citation guidance.
NumPy, SciPy, pybind11, scikit-build-core, threadpoolctl, and h5py are separate
dependencies with their own licenses. This repository's BSD license applies
to its own implementation; dependencies are not relicensed or vendored.
