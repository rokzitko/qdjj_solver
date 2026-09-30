# Bath-accuracy and calculation-time tables

Regenerate with `tools/plot_bath_benchmarks.py`.

Energies and frequencies are in gap units, with $\Delta=1$. Absolute kernel errors and the weak-hybridization Knight-shift slope have units of $1/\Delta$; relative errors are dimensionless. Currents are in units of $2e\Delta/\hbar$. See [the calculation details](https://github.com/rokzitko/qdjj_solver/blob/main/docs/benchmarks/baths/README.md) for reference uncertainties and timing definitions.

## Equal size: eight signed levels per reservoir

| Half-bandwidth | Family | Construction (s) | Max absolute $g$ error | Max relative $g$ error | Max $U=0$ excitation error | Max $U=0$ current error |
|---|---|---|---|---|---|---|
| 10 | cosh-grid | 8.27918e-05 | 0.000535798 | 0.00119072 | 0.00103948 | 2.42738e-06 |
| 10 | surrogate | 0.110544 | 1.39943e-05 | 0.00220185 | 0.00987984 | 8.686e-08 |
| 10 | linear-gl | 4.32921e-05 | 0.308908 | 0.329837 | 0.0947232 | 0.032302 |
| 100 | cosh-grid | 8.02083e-05 | 0.0073311 | 0.0255556 | 0.00451188 | 0.000153847 |
| 100 | surrogate | 0.0329329 | 0.000414273 | 0.0828614 | 0.016038 | 2.95423e-06 |
| 100 | linear-gl | 4.65e-05 | 0.91506 | 0.920923 | 0.534278 | 0.174711 |
| 2000 | cosh-grid | 7.92923e-05 | 0.0292935 | 0.136191 | 0.0138228 | 0.000369034 |
| 2000 | surrogate | 0.0258316 | 0.00104831 | 0.108281 | 0.0173365 | 5.04827e-06 |
| 2000 | linear-gl | 4.4208e-05 | 0.995742 | 0.99606 | 0.934554 | 0.183999 |

The $U=0$ maxima cover all three noninteracting parameter points.

## Surrogate sensitivity at eight levels

| $D$ | Minimum frequency | Maximum fit frequency | Weight exponent | Frequency points | Initial guesses | Random seed | Weak spin slope error | Max $U=0$ excitation error | Max $U=0$ current error | Fit (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.001 | 100 | 0 | 1000 | 4 | 1729 | 1.56164e-05 | 0.016038 | 2.95423e-06 | 0.0329329 |
| 2000 | 0.001 | 100 | 0 | 1000 | 4 | 1729 | 9.12391e-05 | 0.0173365 | 5.04827e-06 | 0.0258316 |
| 100 | 0.001 | 10 | 0 | 1000 | 4 | 1729 | 0.00149348 | 0.0106724 | 1.83193e-05 | 0.0633759 |
| 100 | 0.001 | 1000 | 0 | 1000 | 4 | 1729 | 5.71917e-06 | 0.0162725 | 3.27767e-06 | 0.0328293 |
| 100 | 0.001 | 100 | 0.5 | 1000 | 4 | 1729 | 8.80691e-06 | 0.0192305 | 1.37739e-05 | 0.0225125 |
| 100 | 0.001 | 100 | 1 | 1000 | 4 | 1729 | 1.83527e-06 | 0.0232282 | 6.27888e-05 | 0.0219764 |
| 100 | 1e-06 | 100 | 0 | 1000 | 4 | 1729 | 1.55772e-05 | 0.0159675 | 3.29776e-06 | 0.0367364 |
| 100 | 0.001 | 100 | 0 | 250 | 4 | 1729 | 1.51142e-05 | 0.016048 | 2.97016e-06 | 0.0191299 |
| 100 | 0.001 | 100 | 0 | 4000 | 4 | 1729 | 1.57445e-05 | 0.0160355 | 2.95021e-06 | 0.0941089 |
| 100 | 0.001 | 100 | 0 | 1000 | 1 | 1729 | 1.56164e-05 | 0.016038 | 2.95423e-06 | 0.0083375 |
| 100 | 0.001 | 100 | 0 | 1000 | 8 | 1729 | 1.56164e-05 | 0.016038 | 2.95423e-06 | 0.069919 |
| 100 | 0.001 | 100 | 0 | 1000 | 4 | 1730 | 1.56164e-05 | 0.016038 | 2.95423e-06 | 0.037661 |
| 100 | 0.001 | 100 | 0 | 1000 | 4 | 1731 | 1.56164e-05 | 0.016038 | 2.95423e-06 | 0.0324185 |
| 2000 | 0.001 | 2000 | 0 | 1000 | 4 | 1729 | 1.58559e-05 | 0.0189145 | 9.25093e-06 | 0.0407964 |
| 2000 | 0.001 | 20000 | 0 | 1000 | 4 | 1729 | 1.58348e-05 | 0.0189155 | 9.25501e-06 | 0.0401678 |
| 2000 | 0.001 | 2000 | 1 | 1000 | 4 | 1729 | 0.000458871 | 0.0380469 | 0.00105146 | 0.0281594 |

## Interacting reference checks

| Study | Point | Reference levels | Reference type | Observable | Value | Change with bath size | Change with bond dimension | Difference with another bath family | Max residual |
|---|---|---|---|---|---|---|---|---|---|
| double-dot | strong-p0.5 | 4 | largest finite ED bath | signed_gap | 0.238798 | 0.262454 | — | 0.00174416 | 3.11813e-12 |
| double-dot | strong-p0.9 | 4 | largest finite ED bath | signed_gap | -0.0493728 | 0.437457 | — | 0.00215367 | 1.51848e-12 |
| double-dot | weak-p0.5 | 4 | largest finite ED bath | signed_gap | 0.0365701 | 0.492532 | — | 0.0191943 | 3.4743e-13 |
| double-dot | weak-p0.9 | 4 | largest finite ED bath | signed_gap | 0.00762692 | 0.504814 | — | 0.0333191 | 1.29835e-12 |
| knight | g0.01-p0 | 10 | largest finite ED bath | kappa | 0.00114408 | 2.87124e-06 | — | 3.82419e-07 | 2.79447e-12 |
| knight | g0.01-p0.5 | 10 | largest finite ED bath | kappa | 0.00114488 | 2.88221e-06 | — | 3.8293e-07 | 3.06713e-12 |
| knight | g0.01-p1 | 10 | largest finite ED bath | kappa | 0.00114569 | 2.89318e-06 | — | 3.83442e-07 | 3.1843e-12 |
| knight | g1-p0 | 10 | largest finite ED bath | kappa | 0.134739 | 7.45411e-05 | — | 0.000121148 | 3.19023e-12 |
| knight | g1-p0.5 | 24 | empirical MPS reference | kappa | 0.14329 | 3.57384e-07 | 1.89804e-10 | — | 7.28454e-06 |
| knight | g1-p1 | 10 | largest finite ED bath | kappa | 0.151349 | 8.62421e-05 | — | 2.93828e-05 | 3.01694e-12 |
| knight | g2-p0 | 10 | largest finite ED bath | kappa | 0.297821 | 0.00254408 | — | 0.000566135 | 3.65401e-12 |
| knight | g2-p0.5 | 10 | largest finite ED bath | kappa | 0.325585 | 0.000896261 | — | 0.000305622 | 3.24364e-12 |
| knight | g2-p1 | 10 | largest finite ED bath | kappa | 0.345381 | 0.000149163 | — | 9.3372e-05 | 4.9094e-12 |
| spectrum | g0.5 | 24 | empirical MPS reference | signed_gap | -0.919665 | 5.40287e-05 | 1.33821e-11 | — | 6.14009e-07 |
| spectrum | g10 | 10 | largest finite ED bath | signed_gap | 0.906385 | 0.0038608 | — | 0.00500918 | 8.99682e-13 |
| spectrum | g2 | 10 | largest finite ED bath | signed_gap | -0.292872 | 0.000594796 | — | 1.43198e-06 | 6.26932e-13 |
| spectrum | g20 | 12 | finite MPS reference (bath unresolved) | signed_gap | 0.985462 | — | 3.03348e-08 | — | 2.56153e-05 |
| spectrum | g3 | 10 | largest finite ED bath | signed_gap | 0.0840727 | 0.000353472 | — | 2.61414e-05 | 1.6427e-12 |
| wide-band | p1.50 | 10 | largest finite ED bath | signed_gap | 0.00321637 | 0.0114653 | — | 0.000868312 | 5.72192e-11 |
| wide-band | p1.53 | 10 | largest finite ED bath | signed_gap | -0.00192251 | 0.0114009 | — | 0.000852605 | 4.44787e-11 |
| wide-band | p1.56 | 10 | largest finite ED bath | signed_gap | -0.00713306 | 0.0113366 | — | 0.000836016 | 5.09835e-11 |

## Solver and bath refinements

QP residuals with a finite cutoff are projected residuals; MPS residuals use the declared finite Hamiltonian.

| Study | Point | Bath family | Levels | Method | Representation / ordering | QP cutoff | MPS bond dimension | Value | Finite-problem checks passed | Max residual |
|---|---|---|---|---|---|---|---|---|---|---|
| knight | g1-p0.5 | cosh-grid | 12 | qp | declared-star | 4 | — | 0.143154 | True | 3.66716e-13 |
| knight | g1-p0.5 | cosh-grid | 12 | qp | declared-star | 6 | — | 0.14328 | True | 1.84732e-12 |
| knight | g1-p0.5 | cosh-grid | 12 | dmrg | declared-star | — | 64 | 0.14328 | False | 0.00355257 |
| knight | g1-p0.5 | cosh-grid | 12 | dmrg | declared-star | — | 128 | 0.14328 | False | 0.000301613 |
| knight | g1-p0.5 | cosh-grid | 16 | qp | declared-star | 4 | — | 0.143164 | True | 8.75689e-13 |
| knight | g1-p0.5 | cosh-grid | 16 | qp | declared-star | 6 | — | 0.14329 | True | 3.55295e-12 |
| knight | g1-p0.5 | cosh-grid | 16 | dmrg | declared-star | — | 64 | 0.14329 | False | 0.00440874 |
| knight | g1-p0.5 | cosh-grid | 16 | dmrg | declared-star | — | 128 | 0.14329 | False | 0.000469697 |
| knight | g1-p0.5 | surrogate | 12 | qp | declared-star | 4 | — | 0.143163 | True | 2.84935e-13 |
| knight | g1-p0.5 | surrogate | 12 | qp | declared-star | 6 | — | 0.143289 | True | 8.68386e-13 |
| knight | g1-p0.5 | surrogate | 12 | dmrg | declared-star | — | 64 | 0.143289 | False | 0.00305469 |
| knight | g1-p0.5 | surrogate | 12 | dmrg | declared-star | — | 128 | 0.143289 | False | 0.000240401 |
| knight | g1-p0.5 | surrogate | 16 | qp | declared-star | 4 | — | 0.143164 | True | 7.17745e-13 |
| knight | g1-p0.5 | surrogate | 16 | qp | declared-star | 6 | — | 0.14329 | True | 2.87054e-12 |
| knight | g1-p0.5 | surrogate | 16 | dmrg | declared-star | — | 64 | 0.14329 | False | 0.00384806 |
| knight | g1-p0.5 | surrogate | 16 | dmrg | declared-star | — | 128 | 0.14329 | False | 0.000390287 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | declared-star | — | 64 | -0.919825 | False | 0.00021554 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | declared-star | — | 128 | -0.919825 | False | 1.05195e-05 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | declared-star | — | 64 | -0.919719 | False | 0.000239059 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | declared-star | — | 128 | -0.919719 | False | 1.29789e-05 |
| spectrum | g0.5 | surrogate | 12 | dmrg | declared-star | — | 64 | -0.92277 | False | 0.000178177 |
| spectrum | g0.5 | surrogate | 12 | dmrg | declared-star | — | 128 | -0.92277 | False | 6.97153e-06 |
| spectrum | g20 | cosh-grid | 12 | dmrg | declared-star | — | 64 | 0.985514 | False | 0.0400916 |
| spectrum | g20 | cosh-grid | 12 | dmrg | declared-star | — | 128 | 0.985463 | False | 0.00440331 |
| spectrum | g20 | cosh-grid | 16 | dmrg | declared-star | — | 64 | 0.984477 | False | 0.0437341 |
| spectrum | g20 | cosh-grid | 16 | dmrg | declared-star | — | 128 | 0.984408 | False | 0.00552431 |
| spectrum | g20 | surrogate | 12 | dmrg | declared-star | — | 64 | 0.977726 | False | 0.0279241 |
| spectrum | g20 | surrogate | 12 | dmrg | declared-star | — | 128 | 0.977704 | False | 0.0025754 |
| knight | g1-p0.5 | cosh-grid | 12 | qp | declared-star | 8 | — | 0.14328 | True | 4.00957e-12 |
| knight | g1-p0.5 | cosh-grid | 12 | dmrg | centered | — | 128 | 0.14328 | False | 7.24262e-05 |
| knight | g1-p0.5 | cosh-grid | 12 | dmrg | centered | — | 256 | 0.14328 | True | 2.35708e-06 |
| knight | g1-p0.5 | cosh-grid | 16 | dmrg | centered | — | 128 | 0.14329 | False | 0.00010154 |
| knight | g1-p0.5 | cosh-grid | 16 | dmrg | centered | — | 256 | 0.14329 | True | 5.60378e-06 |
| knight | g1-p0.5 | surrogate | 16 | dmrg | centered | — | 128 | 0.14329 | False | 0.000103489 |
| knight | g1-p0.5 | surrogate | 16 | dmrg | centered | — | 256 | 0.14329 | True | 5.48943e-06 |
| knight | g1-p0.5 | cosh-grid | 24 | dmrg | centered | — | 128 | 0.14329 | False | 0.00011325 |
| knight | g1-p0.5 | cosh-grid | 24 | dmrg | centered | — | 256 | 0.14329 | True | 7.28454e-06 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | electron-chain | — | 128 | -0.919833 | False | 0.0147017 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | electron-chain | — | 256 | -0.919825 | False | 0.000848896 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | electron-chain | — | 128 | -0.919767 | False | 0.0355854 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | electron-chain | — | 256 | -0.919719 | False | 0.00367045 |
| spectrum | g0.5 | cosh-grid | 24 | dmrg | electron-chain | — | 128 | -0.919931 | False | 0.0801274 |
| spectrum | g0.5 | cosh-grid | 24 | dmrg | electron-chain | — | 256 | -0.919671 | False | 0.0124589 |
| spectrum | g20 | cosh-grid | 12 | dmrg | electron-chain | — | 128 | 0.985462 | False | 0.0010224 |
| spectrum | g20 | cosh-grid | 12 | dmrg | electron-chain | — | 256 | 0.985462 | True | 2.56153e-05 |
| spectrum | g20 | cosh-grid | 16 | dmrg | electron-chain | — | 128 | 0.984408 | False | 0.00546239 |
| spectrum | g20 | cosh-grid | 16 | dmrg | electron-chain | — | 256 | 0.984407 | False | 0.000295508 |
| spectrum | g20 | cosh-grid | 24 | dmrg | electron-chain | — | 128 | 0.984239 | False | 0.0222024 |
| spectrum | g20 | cosh-grid | 24 | dmrg | electron-chain | — | 256 | 0.984222 | False | 0.00255306 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | qp-star | — | 128 | -0.919825 | True | 1.05195e-05 |
| spectrum | g0.5 | cosh-grid | 12 | dmrg | qp-star | — | 256 | -0.919825 | True | 2.32326e-07 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | qp-star | — | 128 | -0.919719 | True | 1.29789e-05 |
| spectrum | g0.5 | cosh-grid | 16 | dmrg | qp-star | — | 256 | -0.919719 | True | 4.19977e-07 |
| spectrum | g0.5 | cosh-grid | 24 | dmrg | qp-star | — | 128 | -0.919665 | True | 1.65151e-05 |
| spectrum | g0.5 | cosh-grid | 24 | dmrg | qp-star | — | 256 | -0.919665 | True | 6.14009e-07 |

## Knight shift: g1-p0.5

Reported value: $\kappa$. The point name refers to the parameters in `input.json`.

| Bath family | Levels | Value | Difference from reference | Time with saved bath (s) | Time including bath construction (s) | Peak memory (MiB) |
|---|---|---|---|---|---|---|
| cosh-grid | 2 | 0.176049 | 0.0327585 | 0.00345433 | 0.00351129 | 79.8281 |
| cosh-grid | 4 | 0.141751 | 0.00153929 | 0.00739108 | 0.00745788 | 80.5156 |
| cosh-grid | 6 | 0.141775 | 0.00151497 | 0.039676 | 0.0397486 | 89.8906 |
| cosh-grid | 8 | 0.143448 | 0.000157793 | 0.685841 | 0.685921 | 209.141 |
| cosh-grid | 10 | 0.143356 | 6.65056e-05 | 13.375 | 13.3751 | 1134.66 |
| surrogate | 2 | 0.110995 | 0.032295 | 0.00381796 | 0.0096125 | 81.1719 |
| surrogate | 4 | 0.142057 | 0.00123298 | 0.00722979 | 0.0264978 | 81.4844 |
| surrogate | 6 | 0.143179 | 0.000110621 | 0.0306788 | 0.0541153 | 85.5781 |
| surrogate | 8 | 0.143273 | 1.74371e-05 | 0.504156 | 0.537089 | 207.812 |
| surrogate | 10 | 0.143286 | 4.30582e-06 | 11.2738 | 11.3398 | 1142.61 |
| linear-gl | 2 | 0.0159685 | 0.127321 | 0.00346667 | 0.00350958 | 80.7344 |
| linear-gl | 4 | 0.0293689 | 0.113921 | 0.00695712 | 0.00699821 | 80.4531 |
| linear-gl | 6 | 0.0404888 | 0.102801 | 0.02595 | 0.0259919 | 89.3438 |
| linear-gl | 8 | 0.0498617 | 0.0934283 | 0.372499 | 0.372545 | 201.125 |
| linear-gl | 10 | 0.0578841 | 0.0854058 | 8.34201 | 8.34205 | 1126.38 |

## Parity excitation: g0.5

Reported value: $(E_{\mathrm{odd}}-E_{\mathrm{even}})/\Delta$. The point name refers to the parameters in `input.json`.

| Bath family | Levels | Value | Difference from reference | Time with saved bath (s) | Time including bath construction (s) | Peak memory (MiB) |
|---|---|---|---|---|---|---|
| cosh-grid | 2 | -1.68571 | 0.766049 | 0.00315979 | 0.00321171 | 78.1562 |
| cosh-grid | 4 | -0.970412 | 0.0507474 | 0.00692321 | 0.00699042 | 76.5938 |
| cosh-grid | 6 | -0.914716 | 0.00494925 | 0.0479752 | 0.0480465 | 94.8125 |
| cosh-grid | 8 | -0.915054 | 0.00461077 | 0.873462 | 0.873545 | 304.797 |
| cosh-grid | 10 | -0.918544 | 0.00112126 | 16.8304 | 16.8304 | 1656.66 |
| surrogate | 2 | -1.23284 | 0.313171 | 0.00361808 | 0.00887179 | 81.6562 |
| surrogate | 4 | -1.00661 | 0.0869404 | 0.00725996 | 0.0201129 | 79.7344 |
| surrogate | 6 | -0.952657 | 0.0329917 | 0.0426642 | 0.0656569 | 88.25 |
| surrogate | 8 | -0.933852 | 0.0141871 | 0.827627 | 0.938171 | 330.75 |
| surrogate | 10 | -0.92618 | 0.00651482 | 15.898 | 16.52 | 1674.92 |
| linear-gl | 2 | -4.57623 | 3.65657 | 0.00309596 | 0.00313642 | 78.75 |
| linear-gl | 4 | -2.86794 | 1.94828 | 0.006744 | 0.00678542 | 77.3906 |
| linear-gl | 6 | -2.09974 | 1.18007 | 0.0369796 | 0.0370253 | 82.5312 |
| linear-gl | 8 | -1.6943 | 0.774633 | 0.703072 | 0.703115 | 233 |
| linear-gl | 10 | -1.45512 | 0.535455 | 13.6117 | 13.6118 | 1620.78 |
| surrogate | 1 | -0.803873 | 0.115791 | 0.00252275 | 0.005456 | 77.625 |
| surrogate | 3 | -0.86902 | 0.0506451 | 0.00463775 | 0.0169596 | 78.4844 |
| surrogate | 5 | -0.893919 | 0.0257456 | 0.0141813 | 0.034505 | 77.5781 |
| surrogate | 7 | -0.906114 | 0.0135504 | 0.182353 | 0.241338 | 110.047 |
| surrogate | 9 | -0.912502 | 0.007163 | 3.54815 | 3.7844 | 573.578 |

## Parity excitation: g20

Reported value: $(E_{\mathrm{odd}}-E_{\mathrm{even}})/\Delta$. The point name refers to the parameters in `input.json`.

| Bath family | Levels | Value | Difference from reference | Time with saved bath (s) | Time including bath construction (s) | Peak memory (MiB) |
|---|---|---|---|---|---|---|
| cosh-grid | 2 | 0.912619 | 0.0728439 | 0.00350946 | 0.00356137 | 79.0781 |
| cosh-grid | 4 | 0.963063 | 0.0223997 | 0.00769358 | 0.00776079 | 81.2812 |
| cosh-grid | 6 | 0.977861 | 0.00760195 | 0.0427491 | 0.0428205 | 86.0938 |
| cosh-grid | 8 | 0.983471 | 0.00199126 | 0.74636 | 0.746443 | 267.078 |
| cosh-grid | 10 | 0.985393 | 6.96877e-05 | 14.1943 | 14.1944 | 1607.02 |
| surrogate | 2 | 0.935397 | 0.0500656 | 0.00332754 | 0.00858125 | 77.5781 |
| surrogate | 4 | 0.957073 | 0.028389 | 0.00789388 | 0.0207468 | 81.1719 |
| surrogate | 6 | 0.966488 | 0.0189749 | 0.037478 | 0.0604707 | 94.0781 |
| surrogate | 8 | 0.971896 | 0.0135662 | 0.742352 | 0.852896 | 222.516 |
| surrogate | 10 | 0.975352 | 0.0101106 | 13.4964 | 14.1185 | 1693.56 |
| linear-gl | 2 | 0.423405 | 0.562057 | 0.00322021 | 0.00326067 | 77.7969 |
| linear-gl | 4 | 0.747639 | 0.237824 | 0.00628612 | 0.00632754 | 78.4531 |
| linear-gl | 6 | 0.839942 | 0.14552 | 0.0353139 | 0.0353596 | 94.3906 |
| linear-gl | 8 | 0.881956 | 0.103507 | 0.597749 | 0.597793 | 254.062 |
| linear-gl | 10 | 0.905762 | 0.0797008 | 11.983 | 11.9831 | 1626.62 |
| surrogate | 1 | 2.95535 | 1.96989 | 0.00252892 | 0.00546217 | 79.4062 |
| surrogate | 3 | 1.24347 | 0.258005 | 0.00427425 | 0.0165961 | 78.3906 |
| surrogate | 5 | 1.08948 | 0.10402 | 0.0142547 | 0.0345784 | 82.9844 |
| surrogate | 7 | 1.0374 | 0.0519331 | 0.166906 | 0.22589 | 136.906 |
| surrogate | 9 | 1.01435 | 0.0288901 | 2.95321 | 3.18946 | 641.953 |
