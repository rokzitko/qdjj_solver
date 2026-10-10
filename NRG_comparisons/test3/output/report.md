# Test3: Finite-Field Single-Reservoir Results

## TL;DR

**PILOT ONLY. Production deferred; no continuum extrapolation or uncertainty established.** The full-QP, highest-bond DMRG and highest-retention single-z endpoints are shown first when saved, including the second even root. No empirical NRG uncertainty from Test1 is transferred.

Controls: **passed**. This regenerates saved data, not a new calculation.

## Model

Single superconducting reservoir: gap Delta=1, normal half-bandwidth D=100, U=2, total Gamma=0.4, detuning=0.5, local field B=0.03. Positive B raises spin up. Energies subtract the corresponding disconnected BCS vacuum and retain the centered impurity constant. `even` is parity-even, M=0, not an imposed singlet; the odd states have M=-1/2 and +1/2.

## Energies and Gap

Excitations are relative to each task's lowest saved energy, in gap units. QP and DMRG use the same explicit cosh8 bath; NRG uses Lambda=2, z=1, keep=2000, Nmax=42. The NRG comparison is an unmatched-bath diagnostic, not a reference or z average. Displayed digits are rounded stored values, not established accuracy.

| Endpoint | State | $E/\Delta$ | Excitation / $\Delta$ | Subgap | Degenerate | Residual | Qualification |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cosh8-qp-qfull | odd_down | -1.098284809 | 0 | True | False | 1.532487773e-13 | qualified |
| cosh8-qp-qfull | odd_up | -1.073024061 | 0.02526074803 | True | False | 1.710377731e-13 | qualified |
| cosh8-qp-qfull | even | -1.022332403 | 0.07595240606 | True | False | 1.56818151e-13 | qualified |
| cosh8-qp-qfull | sz+0-r1 | -0.1614129646 | 0.9368718442 | True | False | 1.532232257e-13 | qualified |
| cosh8-dmrg-chi128 | odd_down | -1.098284809 | 0 | True | False | 1.399445507e-06 | not_qualified |
| cosh8-dmrg-chi128 | odd_up | -1.073024061 | 0.02526074803 | True | False | 1.431469849e-06 | not_qualified |
| cosh8-dmrg-chi128 | even | -1.022332403 | 0.07595240606 | True | False | 8.9345016e-06 | not_qualified |
| cosh8-dmrg-chi128 | sz+0-r1 | -0.1614129646 | 0.9368718442 | True | False | 8.204920409e-06 | not_qualified |
| nrg-k2000 | odd_down | -1.103010087 | 0 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | odd_up | -1.077768573 | 0.02524151476 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | even | -1.027512008 | 0.07549807943 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | sz+0-r1 | -0.166892831 | 0.9361172564 | True | False | n/a | unavailable_nrg |

**Additional Subgap States:** the saved pilot exclusion roots include subgap states. `sz+0-r1` labels the second even root, not a discarded observation. This does not establish an additional continuum bound state or spectral completeness.

## Probabilities and Moment

P0/P1/P2 are empty/singly/doubly occupied impurity probabilities. P_up/P_down are singly occupied spin probabilities, not densities. Moments are in hbar units. Each row is state-specific: same-sector degeneracy makes diagonal expectations basis dependent. Qualified full-QP rows concern only this finite bath; DMRG values remain unqualified and NRG residuals are unavailable, not zero.

| Endpoint | State | P0 | P1 | P2 | P_up | P_down | m_dot | m_bath |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cosh8-qp-qfull | odd_down | 0.09068672095 | 0.8580958259 | 0.05121745319 | 0.007347670925 | 0.8507481549 | -0.421700242 | -0.07829975799 |
| cosh8-qp-qfull | odd_up | 0.09237833105 | 0.8558920517 | 0.05172961721 | 0.8482626267 | 0.00762942503 | 0.4203166008 | 0.07968339916 |
| cosh8-qp-qfull | even | 0.539949996 | 0.3692548888 | 0.09079511521 | 0.1812791872 | 0.1879757016 | -0.003348257232 | 0.003348257232 |
| cosh8-qp-qfull | sz+0-r1 | 0.2395031313 | 0.5258267495 | 0.2346701192 | 0.2312706047 | 0.2945561448 | -0.03164277007 | 0.03164277007 |
| cosh8-dmrg-chi128 | odd_down | 0.09068672095 | 0.8580958259 | 0.05121745319 | 0.007347670925 | 0.8507481549 | -0.421700242 | -0.07829975799 |
| cosh8-dmrg-chi128 | odd_up | 0.09237833105 | 0.8558920517 | 0.05172961721 | 0.8482626267 | 0.00762942503 | 0.4203166008 | 0.07968339916 |
| cosh8-dmrg-chi128 | even | 0.539949996 | 0.3692548888 | 0.09079511521 | 0.1812791872 | 0.1879757016 | -0.003348257232 | 0.003348257232 |
| cosh8-dmrg-chi128 | sz+0-r1 | 0.2395031313 | 0.5258267495 | 0.2346701192 | 0.2312706047 | 0.2945561448 | -0.03164277006 | 0.03164277006 |
| nrg-k2000 | odd_down | 0.09105113259 | 0.8577490712 | 0.05119979623 | 0.007474797537 | 0.8502742736 | -0.4213997381 | -0.07860026195 |
| nrg-k2000 | odd_up | 0.09277269267 | 0.8555037052 | 0.05172360215 | 0.8477264924 | 0.007777212738 | 0.4199746398 | 0.08002536015 |
| nrg-k2000 | even | 0.5408054631 | 0.3682368687 | 0.09095766817 | 0.1808651606 | 0.1873717082 | -0.003253273815 | 0.003253273815 |
| nrg-k2000 | sz+0-r1 | 0.2444079866 | 0.5115402783 | 0.2440517351 | 0.2245521348 | 0.2869881435 | -0.03121800435 | 0.03121800435 |

| Endpoint | State | dot_z2 | bath_z2 | Czz | Czz_connected |
| --- | --- | --- | --- | --- | --- |
| cosh8-qp-qfull | odd_down | 0.2145239565 | 0.04282371446 | -0.003673835462 | -0.03669286236 |
| cosh8-qp-qfull | odd_up | 0.2139730129 | 0.0436564121 | -0.003814712515 | -0.03730696799 |
| cosh8-qp-qfull | even | 0.0923137222 | 0.0923137222 | -0.0923137222 | -0.09230251138 |
| cosh8-qp-qfull | sz+0-r1 | 0.1314566874 | 0.1314566874 | -0.1314566874 | -0.1304554225 |
| cosh8-dmrg-chi128 | odd_down | 0.2145239565 | 0.04282371446 | -0.003673835462 | -0.03669286236 |
| cosh8-dmrg-chi128 | odd_up | 0.2139730129 | 0.0436564121 | -0.003814712515 | -0.03730696799 |
| cosh8-dmrg-chi128 | even | 0.0923137222 | 0.0923137222 | -0.0923137222 | -0.09230251138 |
| cosh8-dmrg-chi128 | sz+0-r1 | 0.1314566874 | 0.1314566874 | -0.1314566874 | -0.1304554225 |
| nrg-k2000 | odd_down | 0.2144372678 | 0.04303752974 | -0.003737398768 | -0.03685952856 |
| nrg-k2000 | odd_up | 0.2138759263 | 0.04390128644 | -0.003888606369 | -0.03749722818 |
| nrg-k2000 | even | 0.09205921719 | 0.09205921719 | -0.09205921719 | -0.0920486334 |
| nrg-k2000 | sz+0-r1 | 0.1278850696 | 0.1278850696 | -0.1278850696 | -0.1269105058 |

## Achieved Accuracy

Residual qualification, QP cutoff/MPS bond stability, and bath convergence are separate. Close DMRG energies do not override failed residuals. NRG retention changes at one z and Lambda are not a continuum uncertainty. No backend has established the requested continuum accuracy.

Targets remain absolute 1e-9 plus relative 1e-6 for QP/DMRG; 1e-6 plus 1e-3 for NRG. Normalized changes use the new backend and abs(new). Finite controls instead use their recorded absolute bounds; these targets are not error bars.

## Data and Reproduction

[Results](results.json), [all state observables](observables.csv), and [comparison operands](comparisons.csv) are self-contained. Run `python -B -S NRG_comparisons/test3/generate.py --output-dir PATH` from the repository root. No ignored runs, installed backend, or raw HDF5 file is needed. See the [case README](../README.md) for fresh numerical reproduction.

## Supplemental material

### Refinement Checks

All recorded primary comparisons are retained, including unavailable, unqualified and basis-dependent rows. A small normalized change is not an error bound.

| Kind | Old | New | Compared / Total | Max Normalized Change | Other Statuses |
| --- | --- | --- | --- | --- | --- |
| finite_control | atomic-ed | atomic-nrg | 39/39 | 5.329070518e-09 |  |
| finite_control | atomic-ed | atomic-qp | 39/39 | 2.699229729e-06 |  |
| finite_control | atomic-ed | atomic-dmrg | 39/39 | 2.699229729e-06 |  |
| finite_control | chain-ed | chain-nrg | 39/39 | 1.517399198e-09 |  |
| finite_control | chain-ed | chain-qp | 39/39 | 1.247970726e-06 |  |
| finite_control | chain-ed | chain-dmrg | 39/39 | 1.552438622e-06 |  |
| finite_control | zero-ed | zero-nrg | 39/39 | 1.448263894e-10 |  |
| finite_control | reverse-ed | reverse-nrg | 39/39 | 1.066643992e-09 |  |
| finite_control | chain-nrg | chain-bandwidth | 39/39 | 3.491042806e-09 |  |
| qp_cutoff | cosh8-qp-q2 | cosh8-qp-q4 | 39/39 | 282381.059 |  |
| qp_cutoff | cosh8-qp-q4 | cosh8-qp-q6 | 39/39 | 307.5722846 |  |
| qp_cutoff | cosh8-qp-q6 | cosh8-qp-qfull | 39/39 | 0.02845303952 |  |
| dmrg_bond | cosh8-dmrg-chi32 | cosh8-dmrg-chi64 | 0/39 | n/a | unqualified |
| dmrg_bond | cosh8-dmrg-chi64 | cosh8-dmrg-chi128 | 0/39 | n/a | unqualified |
| nrg_retention | nrg-k500 | nrg-k1000 | 39/39 | 1.708911361 |  |
| nrg_retention | nrg-k1000 | nrg-k2000 | 39/39 | 0.1257285699 |  |
| qp_full_finite_bath | cosh8-qp-q2 | cosh8-qp-qfull | 39/39 | 282747.7659 |  |
| qp_full_finite_bath | cosh8-qp-q4 | cosh8-qp-qfull | 39/39 | 307.6007289 |  |
| dmrg_vs_full_qp | cosh8-qp-qfull | cosh8-dmrg-chi32 | 0/39 | n/a | unqualified |
| dmrg_vs_full_qp | cosh8-qp-qfull | cosh8-dmrg-chi64 | 0/39 | n/a | unqualified |
| dmrg_vs_full_qp | cosh8-qp-qfull | cosh8-dmrg-chi128 | 0/39 | n/a | unqualified |
| unmatched_bath_single_z | cosh8-qp-qfull | nrg-k2000 | 39/39 | 22.33183677 |  |

### Residuals

Every saved root, including exclusions, is retained below and in the CSV. QP residuals refer to the projected cutoff Hamiltonian. DMRG uses the untruncated finite Hamiltonian. An unqualified exclusion root does not automatically invalidate a primary root.

| Task | State | Excitation | Subgap | Degenerate | Residual | Qualification |
| --- | --- | --- | --- | --- | --- | --- |
| atomic-nrg | sz-3-r0 | 109.0758316 | False | False | n/a | unavailable_nrg |
| atomic-nrg | sz-2-r0 | 54.53791579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | sz-2-r1 | 54.53791579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | odd_down | 0 | True | False | n/a | unavailable_nrg |
| atomic-nrg | sz-1-r1 | 55.05291579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | even | 0.515 | True | False | n/a | unavailable_nrg |
| atomic-nrg | sz+0-r1 | 1.515 | False | False | n/a | unavailable_nrg |
| atomic-nrg | odd_up | 0.03 | True | False | n/a | unavailable_nrg |
| atomic-nrg | sz+1-r1 | 55.05291579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | sz+2-r0 | 54.56791579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | sz+2-r1 | 54.56791579 | False | True | n/a | unavailable_nrg |
| atomic-nrg | sz+3-r0 | 109.1058316 | False | False | n/a | unavailable_nrg |
| atomic-ed | sz-3-r0 | 109.0758316 | False | False | 0 | qualified |
| atomic-ed | sz-2-r0 | 54.53791579 | False | True | 5.551115123e-17 | not_qualified |
| atomic-ed | sz-2-r1 | 54.53791579 | False | True | 5.551115123e-17 | not_qualified |
| atomic-ed | odd_down | 0 | True | False | 5.392829661e-14 | qualified |
| atomic-ed | sz-1-r1 | 55.05291579 | False | True | 4.055048304e-14 | not_qualified |
| atomic-ed | even | 0.515 | True | False | 7.666069792e-14 | qualified |
| atomic-ed | sz+0-r1 | 1.515 | False | False | 9.411257969e-14 | qualified |
| atomic-ed | odd_up | 0.03 | True | False | 1.98235407e-14 | qualified |
| atomic-ed | sz+1-r1 | 55.05291579 | False | True | 2.753244355e-14 | not_qualified |
| atomic-ed | sz+2-r0 | 54.56791579 | False | True | 7.106294666e-15 | not_qualified |
| atomic-ed | sz+2-r1 | 54.56791579 | False | True | 7.106294666e-15 | not_qualified |
| atomic-ed | sz+3-r0 | 109.1058316 | False | False | 0 | qualified |
| atomic-qp | sz-3-r0 | 109.0758316 | False | False | 0 | qualified |
| atomic-qp | sz-2-r0 | 54.53791579 | False | True | 0 | not_qualified |
| atomic-qp | sz-2-r1 | 54.53791579 | False | True | 0 | not_qualified |
| atomic-qp | odd_down | 0 | True | False | 0 | qualified |
| atomic-qp | sz-1-r1 | 55.05291579 | False | True | 0 | not_qualified |
| atomic-qp | even | 0.515 | True | False | 0 | qualified |
| atomic-qp | sz+0-r1 | 1.515 | False | False | 0 | qualified |
| atomic-qp | odd_up | 0.03 | True | False | 0 | qualified |
| atomic-qp | sz+1-r1 | 55.05291579 | False | True | 0 | not_qualified |
| atomic-qp | sz+2-r0 | 54.56791579 | False | True | 0 | not_qualified |
| atomic-qp | sz+2-r1 | 54.56791579 | False | True | 0 | not_qualified |
| atomic-qp | sz+3-r0 | 109.1058316 | False | False | 0 | qualified |
| atomic-dmrg | sz-3-r0 | 109.0758316 | False | False | 0 | qualified |
| atomic-dmrg | sz-2-r0 | 54.53791579 | False | True | 0 | not_qualified |
| atomic-dmrg | sz-2-r1 | 54.53791579 | False | True | 0 | not_qualified |
| atomic-dmrg | odd_down | 0 | True | False | 0 | qualified |
| atomic-dmrg | sz-1-r1 | 55.05291579 | False | False | 0 | qualified |
| atomic-dmrg | even | 0.515 | True | False | 2.08307703e-14 | qualified |
| atomic-dmrg | sz+0-r1 | 1.515 | False | False | 0 | qualified |
| atomic-dmrg | odd_up | 0.03 | True | False | 0 | qualified |
| atomic-dmrg | sz+1-r1 | 55.05291579 | False | False | 0 | qualified |
| atomic-dmrg | sz+2-r0 | 54.56791579 | False | True | 0 | not_qualified |
| atomic-dmrg | sz+2-r1 | 54.56791579 | False | True | 0 | not_qualified |
| atomic-dmrg | sz+3-r0 | 109.1058316 | False | False | 0 | qualified |
| chain-nrg | sz-4-r0 | 143.0609232 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz-3-r0 | 74.38027219 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz-3-r1 | 74.38403927 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz-2-r0 | 5.703382615 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz-2-r1 | 71.13075057 | False | False | n/a | unavailable_nrg |
| chain-nrg | odd_down | 2.452652742 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz-1-r1 | 3.673021438 | False | False | n/a | unavailable_nrg |
| chain-nrg | even | 0 | True | False | n/a | unavailable_nrg |
| chain-nrg | sz+0-r1 | 5.598617846 | False | False | n/a | unavailable_nrg |
| chain-nrg | odd_up | 2.471713032 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+1-r1 | 3.688431987 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+2-r0 | 5.733282651 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+2-r1 | 71.14986866 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+3-r0 | 74.41022124 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+3-r1 | 74.41399026 | False | False | n/a | unavailable_nrg |
| chain-nrg | sz+4-r0 | 143.0909232 | False | False | n/a | unavailable_nrg |
| chain-ed | sz-4-r0 | 143.0609232 | False | False | 0 | qualified |
| chain-ed | sz-3-r0 | 74.38027219 | False | False | 3.40698447e-14 | qualified |
| chain-ed | sz-3-r1 | 74.38403927 | False | False | 2.779648841e-14 | qualified |
| chain-ed | sz-2-r0 | 5.703382615 | False | False | 1.688881323e-13 | qualified |
| chain-ed | sz-2-r1 | 71.13075057 | False | False | 7.760298777e-14 | qualified |
| chain-ed | odd_down | 2.452652742 | False | False | 1.387704424e-13 | qualified |
| chain-ed | sz-1-r1 | 3.673021438 | False | False | 1.313588772e-13 | qualified |
| chain-ed | even | 0 | True | False | 1.66044214e-13 | qualified |
| chain-ed | sz+0-r1 | 5.598617846 | False | False | 9.54048687e-14 | qualified |
| chain-ed | odd_up | 2.471713032 | False | False | 1.931752618e-13 | qualified |
| chain-ed | sz+1-r1 | 3.688431987 | False | False | 1.144802322e-13 | qualified |
| chain-ed | sz+2-r0 | 5.733282651 | False | False | 2.667399834e-13 | qualified |
| chain-ed | sz+2-r1 | 71.14986866 | False | False | 8.531920654e-14 | qualified |
| chain-ed | sz+3-r0 | 74.41022124 | False | False | 1.83307454e-14 | qualified |
| chain-ed | sz+3-r1 | 74.41399026 | False | False | 2.62523475e-14 | qualified |
| chain-ed | sz+4-r0 | 143.0909232 | False | False | 0 | qualified |
| chain-qp | sz-4-r0 | 143.0609232 | False | False | 0 | qualified |
| chain-qp | sz-3-r0 | 74.38027219 | False | False | 2.604229767e-14 | qualified |
| chain-qp | sz-3-r1 | 74.38403927 | False | False | 2.864501975e-14 | qualified |
| chain-qp | sz-2-r0 | 5.703382615 | False | False | 1.04466519e-13 | qualified |
| chain-qp | sz-2-r1 | 71.13075057 | False | False | 4.870932156e-14 | qualified |
| chain-qp | odd_down | 2.452652742 | False | False | 2.368903438e-14 | qualified |
| chain-qp | sz-1-r1 | 3.673021438 | False | False | 8.224839948e-14 | qualified |
| chain-qp | even | 0 | True | False | 3.933066852e-14 | qualified |
| chain-qp | sz+0-r1 | 5.598617846 | False | False | 5.954933984e-14 | qualified |
| chain-qp | odd_up | 2.471713032 | False | False | 2.178383866e-14 | qualified |
| chain-qp | sz+1-r1 | 3.688431987 | False | False | 1.120813542e-13 | qualified |
| chain-qp | sz+2-r0 | 5.733282651 | False | False | 1.105965383e-13 | qualified |
| chain-qp | sz+2-r1 | 71.14986866 | False | False | 4.954386667e-14 | qualified |
| chain-qp | sz+3-r0 | 74.41022124 | False | False | 1.881113661e-14 | qualified |
| chain-qp | sz+3-r1 | 74.41399026 | False | False | 3.205070099e-14 | qualified |
| chain-qp | sz+4-r0 | 143.0909232 | False | False | 0 | qualified |
| chain-dmrg | sz-4-r0 | 143.0609232 | False | False | 0 | qualified |
| chain-dmrg | sz-3-r0 | 74.38027219 | False | False | 7.720893507e-14 | qualified |
| chain-dmrg | sz-3-r1 | 74.38403927 | False | False | 9.062023979e-14 | qualified |
| chain-dmrg | sz-2-r0 | 5.703382615 | False | False | 2.494371679e-14 | qualified |
| chain-dmrg | sz-2-r1 | 71.13075057 | False | False | 3.578985578e-11 | qualified |
| chain-dmrg | odd_down | 2.452652742 | False | False | 2.366137768e-14 | qualified |
| chain-dmrg | sz-1-r1 | 3.673021438 | False | False | 2.58917018e-14 | qualified |
| chain-dmrg | even | 0 | True | False | 3.505626159e-14 | qualified |
| chain-dmrg | sz+0-r1 | 5.598617846 | False | False | 6.817960936e-14 | qualified |
| chain-dmrg | odd_up | 2.471713032 | False | False | 3.185070004e-14 | qualified |
| chain-dmrg | sz+1-r1 | 3.688431987 | False | False | 1.993968065e-14 | qualified |
| chain-dmrg | sz+2-r0 | 5.733282651 | False | False | 2.029266223e-14 | qualified |
| chain-dmrg | sz+2-r1 | 71.14986866 | False | False | 7.86388643e-12 | qualified |
| chain-dmrg | sz+3-r0 | 74.41022124 | False | False | 7.948985431e-14 | qualified |
| chain-dmrg | sz+3-r1 | 74.41399026 | False | False | 3.478263957e-14 | qualified |
| chain-dmrg | sz+4-r0 | 143.0909232 | False | False | 0 | qualified |
| zero-nrg | sz-4-r0 | 143.0759027 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz-3-r0 | 74.39522617 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz-3-r1 | 74.39899422 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz-2-r0 | 5.718312091 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz-2-r1 | 71.14029774 | False | False | n/a | unavailable_nrg |
| zero-nrg | odd_down | 2.462171026 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz-1-r1 | 3.680714332 | False | False | n/a | unavailable_nrg |
| zero-nrg | even | 0 | True | False | n/a | unavailable_nrg |
| zero-nrg | sz+0-r1 | 5.598617183 | False | False | n/a | unavailable_nrg |
| zero-nrg | odd_up | 2.462171026 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+1-r1 | 3.680714332 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+2-r0 | 5.718312091 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+2-r1 | 71.14029774 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+3-r0 | 74.39522617 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+3-r1 | 74.39899422 | False | False | n/a | unavailable_nrg |
| zero-nrg | sz+4-r0 | 143.0759027 | False | False | n/a | unavailable_nrg |
| zero-ed | sz-4-r0 | 143.0759027 | False | False | 0 | qualified |
| zero-ed | sz-3-r0 | 74.39522617 | False | False | 3.87901845e-14 | qualified |
| zero-ed | sz-3-r1 | 74.39899422 | False | False | 6.055086205e-14 | qualified |
| zero-ed | sz-2-r0 | 5.718312091 | False | False | 1.479163107e-13 | qualified |
| zero-ed | sz-2-r1 | 71.14029774 | False | False | 6.782837918e-14 | qualified |
| zero-ed | odd_down | 2.462171026 | False | False | 1.900060471e-13 | qualified |
| zero-ed | sz-1-r1 | 3.680714332 | False | False | 1.685115893e-13 | qualified |
| zero-ed | even | 0 | True | False | 1.674922505e-13 | qualified |
| zero-ed | sz+0-r1 | 5.598617183 | False | False | 2.054137449e-13 | qualified |
| zero-ed | odd_up | 2.462171026 | False | False | 1.355676064e-13 | qualified |
| zero-ed | sz+1-r1 | 3.680714332 | False | False | 2.255733686e-13 | qualified |
| zero-ed | sz+2-r0 | 5.718312091 | False | False | 2.196297669e-13 | qualified |
| zero-ed | sz+2-r1 | 71.14029774 | False | False | 4.569718795e-14 | qualified |
| zero-ed | sz+3-r0 | 74.39522617 | False | False | 3.87901845e-14 | qualified |
| zero-ed | sz+3-r1 | 74.39899422 | False | False | 6.055086205e-14 | qualified |
| zero-ed | sz+4-r0 | 143.0759027 | False | False | 0 | qualified |
| reverse-nrg | sz-4-r0 | 143.0909232 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz-3-r0 | 74.41022124 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz-3-r1 | 74.41399026 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz-2-r0 | 5.733282651 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz-2-r1 | 71.14986866 | False | False | n/a | unavailable_nrg |
| reverse-nrg | odd_down | 2.471713032 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz-1-r1 | 3.688431987 | False | False | n/a | unavailable_nrg |
| reverse-nrg | even | 0 | True | False | n/a | unavailable_nrg |
| reverse-nrg | sz+0-r1 | 5.598617846 | False | False | n/a | unavailable_nrg |
| reverse-nrg | odd_up | 2.452652742 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+1-r1 | 3.673021438 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+2-r0 | 5.703382615 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+2-r1 | 71.13075057 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+3-r0 | 74.38027219 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+3-r1 | 74.38403927 | False | False | n/a | unavailable_nrg |
| reverse-nrg | sz+4-r0 | 143.0609232 | False | False | n/a | unavailable_nrg |
| reverse-ed | sz-4-r0 | 143.0909232 | False | False | 0 | qualified |
| reverse-ed | sz-3-r0 | 74.41022124 | False | False | 1.83307454e-14 | qualified |
| reverse-ed | sz-3-r1 | 74.41399026 | False | False | 2.62523475e-14 | qualified |
| reverse-ed | sz-2-r0 | 5.733282651 | False | False | 1.167640303e-13 | qualified |
| reverse-ed | sz-2-r1 | 71.14986866 | False | False | 5.149729521e-14 | qualified |
| reverse-ed | odd_down | 2.471713032 | False | False | 7.547155139e-14 | qualified |
| reverse-ed | sz-1-r1 | 3.688431987 | False | False | 1.174415837e-13 | qualified |
| reverse-ed | even | 0 | True | False | 1.443770743e-13 | qualified |
| reverse-ed | sz+0-r1 | 5.598617846 | False | False | 1.733392696e-13 | qualified |
| reverse-ed | odd_up | 2.452652742 | False | False | 1.572460318e-13 | qualified |
| reverse-ed | sz+1-r1 | 3.673021438 | False | False | 2.304718039e-13 | qualified |
| reverse-ed | sz+2-r0 | 5.703382615 | False | False | 1.527794321e-13 | qualified |
| reverse-ed | sz+2-r1 | 71.13075057 | False | False | 7.495546743e-14 | qualified |
| reverse-ed | sz+3-r0 | 74.38027219 | False | False | 3.40698447e-14 | qualified |
| reverse-ed | sz+3-r1 | 74.38403927 | False | False | 2.779648841e-14 | qualified |
| reverse-ed | sz+4-r0 | 143.0609232 | False | False | 0 | qualified |
| chain-bandwidth | sz-4-r0 | 143.0609232 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz-3-r0 | 74.38027219 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz-3-r1 | 74.38403927 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz-2-r0 | 5.703382615 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz-2-r1 | 71.13075057 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | odd_down | 2.452652742 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz-1-r1 | 3.673021438 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | even | 0 | True | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+0-r1 | 5.598617846 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | odd_up | 2.471713032 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+1-r1 | 3.688431987 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+2-r0 | 5.733282651 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+2-r1 | 71.14986866 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+3-r0 | 74.41022124 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+3-r1 | 74.41399026 | False | False | n/a | unavailable_nrg |
| chain-bandwidth | sz+4-r0 | 143.0909232 | False | False | n/a | unavailable_nrg |
| hf-0.0002--1 | even | 0 | True | False | 1.563892373e-13 | qualified |
| hf-0.0002--1 | sz+0-r1 | 5.598617837 | False | False | 1.632889413e-13 | qualified |
| hf-0.0002--1 | odd_down | 2.452716119 | False | False | 1.974558878e-13 | qualified |
| hf-0.0002--1 | sz-1-r1 | 3.673072642 | False | False | 1.660393359e-13 | qualified |
| hf-0.0002--1 | odd_up | 2.471649341 | False | False | 1.136409756e-13 | qualified |
| hf-0.0002--1 | sz+1-r1 | 3.688380454 | False | False | 1.199239695e-13 | qualified |
| hf-0.0002-+1 | even | 0 | True | False | 1.502587423e-13 | qualified |
| hf-0.0002-+1 | sz+0-r1 | 5.598617855 | False | False | 1.22523481e-13 | qualified |
| hf-0.0002-+1 | odd_down | 2.452589367 | False | False | 1.528517235e-13 | qualified |
| hf-0.0002-+1 | sz-1-r1 | 3.672970236 | False | False | 1.53177505e-13 | qualified |
| hf-0.0002-+1 | odd_up | 2.471776725 | False | False | 2.026643805e-13 | qualified |
| hf-0.0002-+1 | sz+1-r1 | 3.688483521 | False | False | 1.492489456e-13 | qualified |
| hf-0.0001--1 | even | 0 | True | False | 1.192734422e-13 | qualified |
| hf-0.0001--1 | sz+0-r1 | 5.598617842 | False | False | 1.475100514e-13 | qualified |
| hf-0.0001--1 | odd_down | 2.45268443 | False | False | 1.659791354e-13 | qualified |
| hf-0.0001--1 | sz-1-r1 | 3.67304704 | False | False | 1.209535359e-13 | qualified |
| hf-0.0001--1 | odd_up | 2.471681186 | False | False | 1.867992419e-13 | qualified |
| hf-0.0001--1 | sz+1-r1 | 3.68840622 | False | False | 1.639397585e-13 | qualified |
| hf-0.0001-+1 | even | 0 | True | False | 1.719364069e-13 | qualified |
| hf-0.0001-+1 | sz+0-r1 | 5.59861785 | False | False | 1.709501173e-13 | qualified |
| hf-0.0001-+1 | odd_down | 2.452621054 | False | False | 1.448362819e-13 | qualified |
| hf-0.0001-+1 | sz-1-r1 | 3.672995837 | False | False | 1.554990319e-13 | qualified |
| hf-0.0001-+1 | odd_up | 2.471744879 | False | False | 1.802779863e-13 | qualified |
| hf-0.0001-+1 | sz+1-r1 | 3.688457754 | False | False | 1.936252812e-13 | qualified |
| cosh8-qp-q2 | even | 0.09771052695 | True | False | 4.672811628e-14 | qualified |
| cosh8-qp-q2 | sz+0-r1 | 0.9890393406 | True | False | 5.42063434e-14 | qualified |
| cosh8-qp-q2 | odd_down | 0 | True | False | 2.203529077e-14 | qualified |
| cosh8-qp-q2 | sz-1-r1 | 1.442162164 | False | False | 6.301235368e-14 | qualified |
| cosh8-qp-q2 | odd_up | 0.02551508987 | True | False | 1.540602811e-14 | qualified |
| cosh8-qp-q2 | sz+1-r1 | 1.442247298 | False | False | 6.137206227e-14 | qualified |
| cosh8-qp-q2 | sz-2-r0 | 1.172566259 | False | False | 5.560180987e-14 | qualified |
| cosh8-qp-q2 | sz-2-r1 | 1.241684344 | False | False | 6.951108521e-14 | qualified |
| cosh8-qp-q2 | sz+2-r0 | 1.199656665 | False | False | 7.011703897e-14 | qualified |
| cosh8-qp-q2 | sz+2-r1 | 1.269377567 | False | False | 6.414680825e-14 | qualified |
| cosh8-qp-q4 | even | 0.0759744112 | True | False | 7.791801559e-14 | qualified |
| cosh8-qp-q4 | sz+0-r1 | 0.9370214384 | True | False | 1.088945035e-13 | qualified |
| cosh8-qp-q4 | odd_down | 0 | True | False | 1.228901061e-13 | qualified |
| cosh8-qp-q4 | sz-1-r1 | 1.150723794 | False | False | 8.656199039e-14 | qualified |
| cosh8-qp-q4 | odd_up | 0.02526105385 | True | False | 1.45517091e-13 | qualified |
| cosh8-qp-q4 | sz+1-r1 | 1.150831145 | False | False | 1.657034682e-13 | qualified |
| cosh8-qp-q4 | sz-2-r0 | 1.073409951 | False | False | 9.085407518e-14 | qualified |
| cosh8-qp-q4 | sz-2-r1 | 1.170189337 | False | False | 8.520286016e-14 | qualified |
| cosh8-qp-q4 | sz+2-r0 | 1.098783605 | False | False | 8.967039183e-14 | qualified |
| cosh8-qp-q4 | sz+2-r1 | 1.196829198 | False | False | 1.375102474e-13 | qualified |
| cosh8-qp-q6 | even | 0.07595240794 | True | False | 1.086930881e-13 | qualified |
| cosh8-qp-q6 | sz+0-r1 | 0.9368718706 | True | False | 1.838937636e-13 | qualified |
| cosh8-qp-q6 | odd_down | 0 | True | False | 1.575800659e-13 | qualified |
| cosh8-qp-q6 | sz-1-r1 | 1.149436325 | False | False | 1.559238659e-13 | qualified |
| cosh8-qp-q6 | odd_up | 0.02526074806 | True | False | 1.024311885e-13 | qualified |
| cosh8-qp-q6 | sz+1-r1 | 1.149543149 | False | False | 1.211596648e-13 | qualified |
| cosh8-qp-q6 | sz-2-r0 | 1.073163509 | False | False | 1.313478463e-13 | qualified |
| cosh8-qp-q6 | sz-2-r1 | 1.170059989 | False | False | 1.163456545e-13 | qualified |
| cosh8-qp-q6 | sz+2-r0 | 1.098527397 | False | False | 1.623877391e-13 | qualified |
| cosh8-qp-q6 | sz+2-r1 | 1.196695373 | False | False | 1.765714911e-13 | qualified |
| cosh8-qp-qfull | even | 0.07595240606 | True | False | 1.56818151e-13 | qualified |
| cosh8-qp-qfull | sz+0-r1 | 0.9368718442 | True | False | 1.532232257e-13 | qualified |
| cosh8-qp-qfull | odd_down | 0 | True | False | 1.532487773e-13 | qualified |
| cosh8-qp-qfull | sz-1-r1 | 1.149436013 | False | False | 1.861691309e-13 | qualified |
| cosh8-qp-qfull | odd_up | 0.02526074803 | True | False | 1.710377731e-13 | qualified |
| cosh8-qp-qfull | sz+1-r1 | 1.149542837 | False | False | 1.461986167e-13 | qualified |
| cosh8-qp-qfull | sz-2-r0 | 1.073163469 | False | False | 1.930061913e-13 | qualified |
| cosh8-qp-qfull | sz-2-r1 | 1.170059975 | False | False | 1.693609775e-13 | qualified |
| cosh8-qp-qfull | sz+2-r0 | 1.098527356 | False | False | 1.85438942e-13 | qualified |
| cosh8-qp-qfull | sz+2-r1 | 1.196695358 | False | False | 1.323135127e-13 | qualified |
| cosh8-dmrg-chi32 | even | 0.07595511681 | True | False | 0.01572293612 | not_qualified |
| cosh8-dmrg-chi32 | sz+0-r1 | 0.9368753428 | True | False | 0.01747514603 | not_qualified |
| cosh8-dmrg-chi32 | odd_down | 0 | True | False | 0.005275725739 | not_qualified |
| cosh8-dmrg-chi32 | sz-1-r1 | 1.14945225 | False | False | 0.03684019324 | not_qualified |
| cosh8-dmrg-chi32 | odd_up | 0.02526075718 | True | False | 0.005353557884 | not_qualified |
| cosh8-dmrg-chi32 | sz+1-r1 | 1.149558789 | False | False | 0.0365275614 | not_qualified |
| cosh8-dmrg-chi32 | sz-2-r0 | 1.073163506 | False | False | 0.005683626969 | not_qualified |
| cosh8-dmrg-chi32 | sz-2-r1 | 1.17006001 | False | False | 0.005571448601 | not_qualified |
| cosh8-dmrg-chi32 | sz+2-r0 | 1.098527399 | False | False | 0.005742662714 | not_qualified |
| cosh8-dmrg-chi32 | sz+2-r1 | 1.196695401 | False | False | 0.005621501977 | not_qualified |
| cosh8-dmrg-chi64 | even | 0.07595241163 | True | False | 0.0007991169015 | not_qualified |
| cosh8-dmrg-chi64 | sz+0-r1 | 0.9368718507 | True | False | 0.0008923188722 | not_qualified |
| cosh8-dmrg-chi64 | odd_down | 0 | True | False | 0.0002071540791 | not_qualified |
| cosh8-dmrg-chi64 | sz-1-r1 | 1.149436043 | False | False | 0.00179488593 | not_qualified |
| cosh8-dmrg-chi64 | odd_up | 0.02526074805 | True | False | 0.0002121167077 | not_qualified |
| cosh8-dmrg-chi64 | sz+1-r1 | 1.149542866 | False | False | 0.001779299945 | not_qualified |
| cosh8-dmrg-chi64 | sz-2-r0 | 1.073163469 | False | False | 0.00015722853 | not_qualified |
| cosh8-dmrg-chi64 | sz-2-r1 | 1.170059974 | False | False | 0.0001538863838 | not_qualified |
| cosh8-dmrg-chi64 | sz+2-r0 | 1.098527356 | False | False | 0.0001584299517 | not_qualified |
| cosh8-dmrg-chi64 | sz+2-r1 | 1.196695358 | False | False | 0.0001551570004 | not_qualified |
| cosh8-dmrg-chi128 | even | 0.07595240606 | True | False | 8.9345016e-06 | not_qualified |
| cosh8-dmrg-chi128 | sz+0-r1 | 0.9368718442 | True | False | 8.204920409e-06 | not_qualified |
| cosh8-dmrg-chi128 | odd_down | 0 | True | False | 1.399445507e-06 | not_qualified |
| cosh8-dmrg-chi128 | sz-1-r1 | 1.149436013 | False | False | 2.337722992e-05 | not_qualified |
| cosh8-dmrg-chi128 | odd_up | 0.02526074803 | True | False | 1.431469849e-06 | not_qualified |
| cosh8-dmrg-chi128 | sz+1-r1 | 1.149542837 | False | False | 2.310234722e-05 | not_qualified |
| cosh8-dmrg-chi128 | sz-2-r0 | 1.073163469 | False | False | 2.816624363e-07 | not_qualified |
| cosh8-dmrg-chi128 | sz-2-r1 | 1.170059975 | False | False | 2.724070898e-07 | not_qualified |
| cosh8-dmrg-chi128 | sz+2-r0 | 1.098527356 | False | False | 2.855680279e-07 | not_qualified |
| cosh8-dmrg-chi128 | sz+2-r1 | 1.196695358 | False | False | 2.754816511e-07 | not_qualified |
| nrg-k500 | sz-3-r0 | 2.000000004 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz-3-r1 | 2.000000011 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz-2-r0 | 1.000000001 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz-2-r1 | 1.000000004 | False | True | n/a | unavailable_nrg |
| nrg-k500 | odd_down | 0 | True | False | n/a | unavailable_nrg |
| nrg-k500 | sz-1-r1 | 1.075608694 | False | False | n/a | unavailable_nrg |
| nrg-k500 | even | 0.07560869429 | True | False | n/a | unavailable_nrg |
| nrg-k500 | sz+0-r1 | 0.9366538413 | True | False | n/a | unavailable_nrg |
| nrg-k500 | odd_up | 0.02527225609 | True | False | n/a | unavailable_nrg |
| nrg-k500 | sz+1-r1 | 1.075608694 | False | False | n/a | unavailable_nrg |
| nrg-k500 | sz+2-r0 | 1.025272257 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz+2-r1 | 1.02527226 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz+3-r0 | 2.02527226 | False | True | n/a | unavailable_nrg |
| nrg-k500 | sz+3-r1 | 2.025272267 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz-3-r0 | 2.000000004 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz-3-r1 | 2.000000011 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz-2-r0 | 1.000000001 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz-2-r1 | 1.000000004 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | odd_down | 0 | True | False | n/a | unavailable_nrg |
| nrg-k1000 | sz-1-r1 | 1.075488462 | False | False | n/a | unavailable_nrg |
| nrg-k1000 | even | 0.07548846144 | True | False | n/a | unavailable_nrg |
| nrg-k1000 | sz+0-r1 | 0.9362108465 | True | False | n/a | unavailable_nrg |
| nrg-k1000 | odd_up | 0.02524345982 | True | False | n/a | unavailable_nrg |
| nrg-k1000 | sz+1-r1 | 1.075488462 | False | False | n/a | unavailable_nrg |
| nrg-k1000 | sz+2-r0 | 1.025243461 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz+2-r1 | 1.025243464 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz+3-r0 | 2.025243464 | False | True | n/a | unavailable_nrg |
| nrg-k1000 | sz+3-r1 | 2.025243471 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-4-r0 | 3.000000008 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-4-r1 | 3.000000013 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-3-r0 | 2.000000003 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-3-r1 | 2.000000005 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-2-r0 | 1.000000001 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz-2-r1 | 1.000000002 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | odd_down | 0 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | sz-1-r1 | 1.07549808 | False | False | n/a | unavailable_nrg |
| nrg-k2000 | even | 0.07549807943 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | sz+0-r1 | 0.9361172564 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | odd_up | 0.02524151476 | True | False | n/a | unavailable_nrg |
| nrg-k2000 | sz+1-r1 | 1.07549808 | False | False | n/a | unavailable_nrg |
| nrg-k2000 | sz+2-r0 | 1.025241516 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz+2-r1 | 1.025241519 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz+3-r0 | 2.025241519 | False | True | n/a | unavailable_nrg |
| nrg-k2000 | sz+3-r1 | 2.025241525 | False | True | n/a | unavailable_nrg |

### Controls


| Recorded Control | Maximum Error | Tolerance | Passed |
| --- | --- | --- | --- |
| atomic-ed declared root inventory | n/a | n/a | True |
| chain-ed declared root inventory | n/a | n/a | True |
| zero-ed declared root inventory | n/a | n/a | True |
| reverse-ed declared root inventory | n/a | n/a | True |
| atomic-nrg covers all atomic-ed control roots | n/a | n/a | True |
| atomic-nrg vs atomic-ed | 1.532107774e-13 | 1e-08 | True |
| atomic-qp covers all atomic-ed control roots | n/a | n/a | True |
| atomic-qp vs atomic-ed | 8.215650382e-14 | 1e-08 | True |
| atomic-dmrg covers all atomic-ed control roots | n/a | n/a | True |
| atomic-dmrg vs atomic-ed | 8.215650382e-14 | 1e-08 | True |
| chain-nrg covers all chain-ed control roots | n/a | n/a | True |
| chain-nrg vs chain-ed | 1.065814104e-13 | 1e-08 | True |
| chain-qp covers all chain-ed control roots | n/a | n/a | True |
| chain-qp vs chain-ed | 1.07913678e-13 | 1e-08 | True |
| chain-dmrg covers all chain-ed control roots | n/a | n/a | True |
| chain-dmrg vs chain-ed | 7.838174554e-14 | 1e-08 | True |
| zero-nrg covers all zero-ed control roots | n/a | n/a | True |
| zero-nrg vs zero-ed | 1.779687508e-13 | 1e-08 | True |
| reverse-nrg covers all reverse-ed control roots | n/a | n/a | True |
| reverse-nrg vs reverse-ed | 7.815970093e-14 | 1e-08 | True |
| chain-bandwidth covers all chain-nrg control roots | n/a | n/a | True |
| chain-bandwidth vs chain-nrg | 1.136868377e-13 | 1e-08 | True |
| atomic energies and Zeeman sign | 7.105427358e-14 | 1e-09 | True |
| NRG spin reversal | 8.526512829e-14 | 1e-08 | True |
| Hellmann-Feynman h=0.0002 | 3.895975209e-10 | 1e-07 | True |
| Hellmann-Feynman h=0.0001 | 2.286949519e-10 | 1e-07 | True |
| atomic-ed primary residual/uniqueness qualification | n/a | n/a | True |
| atomic-ed all control-root residuals and orthogonality | n/a | n/a | True |
| atomic-qp primary residual/uniqueness qualification | n/a | n/a | True |
| atomic-qp all control-root residuals and orthogonality | n/a | n/a | True |
| atomic-dmrg primary residual/uniqueness qualification | n/a | n/a | True |
| atomic-dmrg all control-root residuals and orthogonality | n/a | n/a | True |
| chain-ed primary residual/uniqueness qualification | n/a | n/a | True |
| chain-ed all control-root residuals and orthogonality | n/a | n/a | True |
| chain-qp primary residual/uniqueness qualification | n/a | n/a | True |
| chain-qp all control-root residuals and orthogonality | n/a | n/a | True |
| chain-dmrg primary residual/uniqueness qualification | n/a | n/a | True |
| chain-dmrg all control-root residuals and orthogonality | n/a | n/a | True |
| zero-ed primary residual/uniqueness qualification | n/a | n/a | True |
| zero-ed all control-root residuals and orthogonality | n/a | n/a | True |
| reverse-ed primary residual/uniqueness qualification | n/a | n/a | True |
| reverse-ed all control-root residuals and orthogonality | n/a | n/a | True |

### Task Status and Costs

Elapsed task times and memory are measurements, not scaling estimates.

| Task | Stage | Status | Primary Qualified | Seconds | Peak GiB | Issue |
| --- | --- | --- | --- | --- | --- | --- |
| atomic-nrg | controls | completed | n/a | 3.603640068 | 0.2932319641 |  |
| atomic-ed | controls | completed | True | 0.2852879986 | 0.01103591919 |  |
| atomic-qp | controls | completed | True | 0.5604212992 | 0.06910324097 |  |
| atomic-dmrg | controls | completed | True | 5.248580784 | 0.1193008423 |  |
| chain-nrg | controls | completed | n/a | 3.047029808 | 0.2350769043 |  |
| chain-ed | controls | completed | True | 0.2854547203 | 0.01103591919 |  |
| chain-qp | controls | completed | True | 0.5632813051 | 0.06923294067 |  |
| chain-dmrg | controls | completed | True | 13.78727802 | 0.1187667847 |  |
| zero-nrg | controls | completed | n/a | 3.049137082 | 0.2353401184 |  |
| zero-ed | controls | completed | True | 0.2857418135 | 0.0111579895 |  |
| reverse-nrg | controls | completed | n/a | 2.770535793 | 0.2357597351 |  |
| reverse-ed | controls | completed | True | 0.2861562148 | 0.01091003418 |  |
| chain-bandwidth | controls | completed | n/a | 2.769532751 | 0.2361755371 |  |
| hf-0.0002--1 | controls | completed | True | 0.2854730524 | 0.01103591919 |  |
| hf-0.0002-+1 | controls | completed | True | 0.2853749394 | 0.01121902466 |  |
| hf-0.0001--1 | controls | completed | True | 0.2826122679 | 0.01103591919 |  |
| hf-0.0001-+1 | controls | completed | True | 0.2825098597 | 0.01098632812 |  |
| cosh8-qp-q2 | pilot | completed | True | 0.5616350509 | 0.06918334961 |  |
| cosh8-qp-q4 | pilot | completed | True | 0.8383898437 | 0.06899261475 |  |
| cosh8-qp-q6 | pilot | completed | True | 1.38964957 | 0.08243179321 |  |
| cosh8-qp-qfull | pilot | completed | True | 5.265043251 | 0.1244239807 |  |
| cosh8-dmrg-chi32 | pilot | completed | False | 456.760884 | 0.1302986145 |  |
| cosh8-dmrg-chi64 | pilot | completed | False | 428.7498195 | 0.147277832 |  |
| cosh8-dmrg-chi128 | pilot | completed | False | 263.8152126 | 0.1712341309 |  |
| nrg-k500 | pilot | completed | n/a | 3.60482017 | 0.3058166504 |  |
| nrg-k1000 | pilot | completed | n/a | 5.25191392 | 0.3593025208 |  |
| nrg-k2000 | pilot | completed | n/a | 14.66456488 | 0.5226325989 |  |

### Limitations

- Three primary branches are selected, not proof of spectral completeness. expected_three_subgap_states describes returned roots only.
- Excitations use each task's lowest saved energy; untested sectors and roots can change that minimum.
- Same-sector degeneracies can make diagonal expectations basis dependent; cross-sector energy ties do not.
- Finite-problem residual/uniqueness qualification is separate from QP cutoff, MPS bond and bath convergence.
- NRG residuals are unavailable. Retention changes at fixed Lambda, z and chain length are not continuum uncertainty.
- Targets normalize differences by absolute + relative*abs(new), using the new backend (ED uses QP targets).
- Finite controls use a separate absolute 1e-8 comparison; recorded analytic and derivative checks have their own tolerances.
- Single-z NRG versus cosh8 full QP is an unmatched-bath diagnostic, not a reference, z average or uncertainty estimate.
- No uncertainty or convergence claim is inherited from test1, test2, or another campaign. Production requires separate approval.

### Provenance

The public projection keeps explicit baths, chain coefficients, numerical inputs, all saved state values, control checks, status and qualifications. Original receipt/source hashes identify archived originals, not the projected objects or current source files. Public request/result hashes separately identify sorted-key JSON of the projected nested objects. Raw artifact hashes are copied, not revalidated here. Full NRG spectra remain in the original receipts; compact sector summaries retain omitted roots and gap-edge evidence. Private execution identities are not publication data.
