# Test4 Working Results

Coarse working evidence, not a continuum reference or error bar. Production approval, continuum accuracy and spectral completeness remain **false**.

The separate mixed-retention Lambda8 Nz8 mean uses seven keep4000 points and keep5000 only at z=0.125. Both profiles use Nmax=16, keepenergy=16 and keepmin=400, with identical twists in the two physical leads. The original uniform keep4000 mean is unchanged. Finite cosh12 QP q6 is a bath- and cutoff-unconverged comparator, not truth.

Energies/excitations are in gap units; current is $I/(2e\Delta/\hbar)$. Gamma=0.4 is total hybridization, D=100 is half-bandwidth, U=2, nu=0.75, B=0.03, phi=3*pi/5, WN=WS=0. The paper energy offset is 0.0625.

Requested comparison tolerances remain 1e-9 + 1e-6*abs(new) for QP/DMRG and 1e-6 + 0.001*abs(new) for NRG. The reported discrepancies and 1%/5% diagnostic tiers are not independently established continuum error bars.

| State (parity, twice Sz, root) | Quantity | Uniform NRG | Mixed NRG | QP cosh12 q6 | Discrepancy % |
| --- | --- | --- | --- | --- | --- |
| (1, -1, 0) | energy | -1.098614553 | -1.105084686 | -1.111336132 | 0.5656984386 |
| (1, -1, 0) | P0 | 0.09614116455 | 0.0963320783 | 0.09607791288 | 0.2638429745 |
| (1, -1, 0) | P1 | 0.8497023652 | 0.8494535798 | 0.8498900165 | 0.05137852111 |
| (1, -1, 0) | P2 | 0.05415647027 | 0.0542143419 | 0.05403207064 | 0.3362048853 |
| (1, -1, 0) | P_up | 0.01051772479 | 0.01050690923 | 0.01044953723 | 0.5460406822 |
| (1, -1, 0) | P_down | 0.8391846404 | 0.8389466706 | 0.8394404793 | 0.05886055728 |
| (1, -1, 0) | m_dot | -0.4143334578 | -0.4142198807 | -0.414495471 | 0.06653237931 |
| (1, -1, 0) | m_bath | -0.08566654219 | -0.08578011933 | -0.08550452898 | 0.3212753091 |
| (1, -1, 0) | dot_z2 | 0.2124255913 | 0.2123633949 | 0.2124725041 | 0.05137852127 |
| (1, -1, 0) | bath_z2 | 0.04809213349 | 0.04814351428 | 0.04797703311 | 0.3458018598 |
| (1, -1, 0) | Czz | -0.005258862393 | -0.005253454613 | -0.005224768617 | 0.5460406124 |
| (1, -1, 0) | Czz_connected | -0.04075268317 | -0.04078459205 | -0.04066600863 | 0.2907554376 |
| (1, -1, 0) | energy_paper | -1.036114553 | -1.042584686 | -1.048836132 | 0.5996104583 |
| (1, -1, 0) | current | -0.00933012384 | -0.00941713864 | -0.009397873601 | 0.2045742334 |
| (1, -1, 0) | excitation | 0 | 0 | 0 | n/a |
| (1, 1, 0) | energy | -1.073621681 | -1.079997766 | -1.086519303 | 0.6038472251 |
| (1, 1, 0) | P0 | 0.09808789772 | 0.09832562601 | 0.09801996005 | 0.3108711082 |
| (1, 1, 0) | P1 | 0.847125902 | 0.8468236079 | 0.8473216404 | 0.05881183053 |
| (1, 1, 0) | P2 | 0.05478620031 | 0.05485076608 | 0.05465839958 | 0.35070887 |
| (1, 1, 0) | P_up | 0.8361109605 | 0.8358178354 | 0.8363803906 | 0.06730595155 |
| (1, 1, 0) | P_down | 0.01101494143 | 0.01100577247 | 0.01094124979 | 0.5862621849 |
| (1, 1, 0) | m_dot | 0.4125480096 | 0.4124060315 | 0.4127195704 | 0.07602675296 |
| (1, 1, 0) | m_bath | 0.08745199044 | 0.08759396852 | 0.0872804296 | 0.3579457831 |
| (1, 1, 0) | dot_z2 | 0.2117814755 | 0.211705902 | 0.2118304101 | 0.05881183072 |
| (1, 1, 0) | bath_z2 | 0.04923346594 | 0.0492998705 | 0.0491108397 | 0.3834306129 |
| (1, 1, 0) | Czz | -0.005507470714 | -0.005502886237 | -0.0054706249 | 0.5862621175 |
| (1, 1, 0) | Czz_connected | -0.04158485117 | -0.04162637128 | -0.04149296631 | 0.3204818597 |
| (1, 1, 0) | energy_paper | -1.011121681 | -1.017497766 | -1.024019303 | 0.6409386594 |
| (1, 1, 0) | current | -0.009627038646 | -0.009716108073 | -0.009689263144 | 0.2762930182 |
| (1, 1, 0) | excitation | 0.02499287141 | 0.0250869197 | 0.02481682997 | 1.076615752 |
| (1, 1, 0) | odd_splitting | 0.02499287141 | 0.0250869197 | 0.02481682997 | 1.076615752 |
| (0, 0, 0) | energy | -0.9374953052 | -0.9470135954 | -0.9532064665 | 0.6539368759 |
| (0, 0, 0) | P0 | 0.5958128612 | 0.5984722699 | 0.5994268251 | 0.1594986482 |
| (0, 0, 0) | P1 | 0.3476921335 | 0.3454979501 | 0.3446018322 | 0.2593699703 |
| (0, 0, 0) | P2 | 0.05649500523 | 0.05602978 | 0.05597134275 | 0.1042967687 |
| (0, 0, 0) | P_up | 0.1705340403 | 0.1694751773 | 0.1690649454 | 0.2420601362 |
| (0, 0, 0) | P_down | 0.1771580932 | 0.1760227728 | 0.1755368867 | 0.276035923 |
| (0, 0, 0) | m_dot | -0.003312026466 | -0.003273797788 | -0.003235970668 | 1.155450729 |
| (0, 0, 0) | m_bath | 0.003312026466 | 0.003273797788 | 0.003235970668 | 1.155450729 |
| (0, 0, 0) | dot_z2 | 0.08692303339 | 0.08637448752 | 0.08615045804 | 0.2593699693 |
| (0, 0, 0) | bath_z2 | 0.08692303339 | 0.08637448752 | 0.08615045804 | 0.2593699693 |
| (0, 0, 0) | Czz | -0.08692303339 | -0.08637448752 | -0.08615045804 | 0.2593699693 |
| (0, 0, 0) | Czz_connected | -0.08691198776 | -0.08636373051 | -0.08613998653 | 0.2590716864 |
| (0, 0, 0) | energy_paper | -0.8749953052 | -0.8845135954 | -0.8907064665 | 0.7001442547 |
| (0, 0, 0) | current | 0.06013131459 | 0.0596542933 | 0.05956223761 | 0.1543152759 |
| (0, 0, 0) | excitation | 0.1611192476 | 0.1580710903 | 0.158129666 | 0.03705654964 |
| (0, 0, 1) | energy | -0.330967703 | -0.3394229158 | -0.345980403 | 1.931951817 |
| (0, 0, 1) | P0 | 0.1422378462 | 0.1419840792 | 0.1415345601 | 0.3165981919 |
| (0, 0, 1) | P1 | 0.5759269791 | 0.5736571388 | 0.5734849496 | 0.03001603857 |
| (0, 0, 1) | P2 | 0.2818351747 | 0.284358782 | 0.2849804902 | 0.2186351239 |
| (0, 0, 1) | P_up | 0.2707777986 | 0.2697998708 | 0.2695419851 | 0.09558407086 |
| (0, 0, 1) | P_down | 0.3051491805 | 0.3038572679 | 0.3039429645 | 0.02820289681 |
| (0, 0, 1) | m_dot | -0.01718569095 | -0.01702869854 | -0.01720048967 | 1.008832973 |
| (0, 0, 1) | m_bath | 0.01718569095 | 0.01702869854 | 0.01720048967 | 1.008832973 |
| (0, 0, 1) | dot_z2 | 0.1439817448 | 0.1434142847 | 0.1433712374 | 0.03001603736 |
| (0, 0, 1) | bath_z2 | 0.1439817448 | 0.1434142847 | 0.1433712374 | 0.03001603736 |
| (0, 0, 1) | Czz | -0.1439817448 | -0.1434142847 | -0.1433712374 | 0.03001603736 |
| (0, 0, 1) | Czz_connected | -0.1436836664 | -0.1431221227 | -0.1430753806 | 0.03265890107 |
| (0, 0, 1) | energy_paper | -0.268467703 | -0.2769229158 | -0.283480403 | 2.367982862 |
| (0, 0, 1) | current | -0.11982534 | -0.1195411001 | -0.119036787 | 0.4218742178 |
| (0, 0, 1) | excitation | 0.7676468497 | 0.7656617699 | 0.7653557295 | 0.03997070094 |

Discrepancy = 100*abs(NRG-QP)/abs(NRG); near zero (<=1e-6), only absolute differences are meaningful. These are empirical deltas, not statistical errors. Absolute-energy percentages depend on the stated energy zero. Each excitation and connected Czz is formed at each z before averaging.

## Supplemental material

### Controls

The 37 original controls retain all returned roots, including above-gap and degenerate rows. Atomic, unit, representation, phase-Hellmann-Feynman and reversal checks below are scalar-replayed. Solver norm/Gram/sweep/root checks are recorded evidence, not newly checked raw MPS/HDF5.

| Check | Maximum Error | Tolerance | Passed |
| --- | --- | --- | --- |
| atomic-ed inventory | 0 | 0.5 | True |
| atomic-ed recorded finite checks | 0 | 0.5 | True |
| atomic-qp inventory | 0 | 0.5 | True |
| atomic-qp recorded finite checks | 0 | 0.5 | True |
| bridge-ed inventory | 0 | 0.5 | True |
| bridge-ed recorded finite checks | 0 | 0.5 | True |
| bridge-qp inventory | 0 | 0.5 | True |
| bridge-qp recorded finite checks | 0 | 0.5 | True |
| direct-ed inventory | 0 | 0.5 | True |
| direct-ed recorded finite checks | 0 | 0.5 | True |
| free-ed inventory | 0 | 0.5 | True |
| free-ed recorded finite checks | 0 | 0.5 | True |
| free-qp inventory | 0 | 0.5 | True |
| free-qp recorded finite checks | 0 | 0.5 | True |
| hf-bridge-0.0001-+1 inventory | 0 | 0.5 | True |
| hf-bridge-0.0001-+1 recorded finite checks | 0 | 0.5 | True |
| hf-bridge-0.0001--1 inventory | 0 | 0.5 | True |
| hf-bridge-0.0001--1 recorded finite checks | 0 | 0.5 | True |
| hf-bridge-0.0002-+1 inventory | 0 | 0.5 | True |
| hf-bridge-0.0002-+1 recorded finite checks | 0 | 0.5 | True |
| hf-bridge-0.0002--1 inventory | 0 | 0.5 | True |
| hf-bridge-0.0002--1 recorded finite checks | 0 | 0.5 | True |
| hf-direct-0.0001-+1 inventory | 0 | 0.5 | True |
| hf-direct-0.0001-+1 recorded finite checks | 0 | 0.5 | True |
| hf-direct-0.0001--1 inventory | 0 | 0.5 | True |
| hf-direct-0.0001--1 recorded finite checks | 0 | 0.5 | True |
| hf-direct-0.0002-+1 inventory | 0 | 0.5 | True |
| hf-direct-0.0002-+1 recorded finite checks | 0 | 0.5 | True |
| hf-direct-0.0002--1 inventory | 0 | 0.5 | True |
| hf-direct-0.0002--1 recorded finite checks | 0 | 0.5 | True |
| hf-junction-0.0001-+1 inventory | 0 | 0.5 | True |
| hf-junction-0.0001-+1 recorded finite checks | 0 | 0.5 | True |
| hf-junction-0.0001--1 inventory | 0 | 0.5 | True |
| hf-junction-0.0001--1 recorded finite checks | 0 | 0.5 | True |
| hf-junction-0.0002-+1 inventory | 0 | 0.5 | True |
| hf-junction-0.0002-+1 recorded finite checks | 0 | 0.5 | True |
| hf-junction-0.0002--1 inventory | 0 | 0.5 | True |
| hf-junction-0.0002--1 recorded finite checks | 0 | 0.5 | True |
| junction-dmrg inventory | 0 | 0.5 | True |
| junction-dmrg recorded finite checks | 0 | 0.5 | True |
| junction-ed inventory | 0 | 0.5 | True |
| junction-ed recorded finite checks | 0 | 0.5 | True |
| junction-qp inventory | 0 | 0.5 | True |
| junction-qp recorded finite checks | 0 | 0.5 | True |
| normal-ed inventory | 0 | 0.5 | True |
| normal-ed recorded finite checks | 0 | 0.5 | True |
| phase-reverse-ed inventory | 0 | 0.5 | True |
| phase-reverse-ed recorded finite checks | 0 | 0.5 | True |
| reverse-ed inventory | 0 | 0.5 | True |
| reverse-ed recorded finite checks | 0 | 0.5 | True |
| spin-ed inventory | 0 | 0.5 | True |
| spin-ed recorded finite checks | 0 | 0.5 | True |
| junction-bandwidth covers junction-nrg | 0 | 0.5 | True |
| junction-bandwidth vs junction-nrg | 1.421085472e-13 | 1e-08 | True |
| junction-sds covers junction-nrg | 0 | 0.5 | True |
| junction-sds vs junction-nrg | 1.98951966e-13 | 1e-08 | True |
| atomic-nrg covers atomic-ed | 0 | 0.5 | True |
| atomic-nrg vs atomic-ed | 1.736388811e-13 | 1e-08 | True |
| atomic-qp covers atomic-ed | 0 | 0.5 | True |
| atomic-qp vs atomic-ed | 1.776356839e-13 | 1e-08 | True |
| bridge-nrg covers bridge-ed | 0 | 0.5 | True |
| bridge-nrg vs bridge-ed | 5.359931249e-13 | 1e-08 | True |
| bridge-qp covers bridge-ed | 0 | 0.5 | True |
| bridge-qp vs bridge-ed | 4.083729882e-13 | 1e-08 | True |
| direct-nrg covers direct-ed | 0 | 0.5 | True |
| direct-nrg vs direct-ed | 3.126388037e-13 | 1e-08 | True |
| free-nrg covers free-ed | 0 | 0.5 | True |
| free-nrg vs free-ed | 9.865441797e-13 | 1e-08 | True |
| free-qp covers free-ed | 0 | 0.5 | True |
| free-qp vs free-ed | 1.08446585e-12 | 1e-08 | True |
| junction-dmrg covers junction-ed | 0 | 0.5 | True |
| junction-dmrg vs junction-ed | 1.207922651e-13 | 1e-08 | True |
| junction-nrg covers junction-ed | 0 | 0.5 | True |
| junction-nrg vs junction-ed | 1.98951966e-13 | 1e-08 | True |
| junction-qp covers junction-ed | 0 | 0.5 | True |
| junction-qp vs junction-ed | 2.629008122e-13 | 1e-08 | True |
| normal-nrg covers normal-ed | 0 | 0.5 | True |
| normal-nrg vs normal-ed | 2.682298827e-13 | 1e-08 | True |
| phase-reverse-nrg covers phase-reverse-ed | 0 | 0.5 | True |
| phase-reverse-nrg vs phase-reverse-ed | 5.359931249e-13 | 1e-08 | True |
| reverse-nrg covers reverse-ed | 0 | 0.5 | True |
| reverse-nrg vs reverse-ed | 1.88293825e-13 | 1e-08 | True |
| spin-nrg covers spin-ed | 0 | 0.5 | True |
| spin-nrg vs spin-ed | 1.867395127e-13 | 1e-08 | True |
| atomic energies, paper offset, Zeeman sign, current | 1.705302566e-13 | 1e-09 | True |
| junction phase HF h=0.0002 | 3.316852323e-10 | 1e-07 | True |
| junction phase HF h=0.0001 | 4.31050193e-10 | 1e-07 | True |
| bridge phase HF h=0.0002 | 2.402721676e-10 | 1e-07 | True |
| bridge phase HF h=0.0001 | 8.477463965e-10 | 1e-07 | True |
| direct phase HF h=0.0002 | 5.564398131e-10 | 1e-07 | True |
| direct phase HF h=0.0001 | 8.557589862e-10 | 1e-07 | True |
| Gamma=0 direct bridge current nonzero | 0 | 0.5 | True |
| reverse symmetry | 1.136868377e-13 | 1e-08 | True |
| phase-reverse symmetry | 0 | 1e-08 | True |

Lambda8 untruncated control: 1024 energies per backend (NRG/ED/full QP); 576 multiplets, 256 unique. Maximum scalar-replayed error: 4.206579529e-12. Degenerate individual expectations are excluded; ED/QP multiplet means are compared. Both actual Wilson leads and unnormalized star weights are retained.

### Convergence

Requested targets use the new backend: QP/DMRG 1e-9 + 1e-6*abs(new), NRG 1e-6 + 0.001*abs(new). ED controls retain their separate absolute tolerance. The additional 1%/5% diagnostics use 1e-6 + relative*abs(new); they do not replace a finite-backend target. Passing a step is not a convergence certificate.

- endpoint_cap_4000_to_5000: 1 / 18 / 38 of 60 scalars (requested / 1% / 5%).
- endpoint_mean_impact: 15 / 53 / 60 of 60 scalars (requested / 1% / 5%).
- finite_bath_QP_diagnostic: 18 / 53 / 60 of 60 scalars (requested / 1% / 5%).

QP residuals qualify only the projected finite problem, not cutoff or bath convergence. DMRG qualification requires recorded finite-converged, norm, Gram, sweep and same-call root checks. Single-sector calls have no global excitations; compatible assemblies must share physics, bath, full Hamiltonian fingerprint and coordinates, with no even-root splicing.


| Convergence Diagnostic | Scalar Pairs | Qualified | Requested | 1% | 5% |
| --- | --- | --- | --- | --- | --- |
| cross_lambda_discretization_and_finite_bath | 120 | 60 | 9 | 48 | 58 |
| dmrg_bond | 532 | 0 | 0 | 0 | 0 |
| dmrg_state_vs_qp | 826 | 84 | 84 | 84 | 84 |
| dmrg_vs_full_qp | 120 | 120 | 120 | 120 | 120 |
| dmrg_vs_qp | 240 | 15 | 15 | 15 | 15 |
| finite_bath_size | 60 | 45 | 1 | 1 | 7 |
| nrg_Nz | 480 | 420 | 135 | 325 | 410 |
| nrg_keep | 1560 | 1500 | 761 | 1155 | 1397 |
| nrg_keepenergy | 300 | 300 | 300 | 300 | 300 |
| nrg_keepmin | 60 | 60 | 60 | 60 | 60 |
| nrg_mean_cap | 360 | 300 | 152 | 269 | 297 |
| nrg_nmax | 360 | 300 | 300 | 300 | 300 |
| qp_cutoff | 600 | 600 | 246 | 547 | 580 |

Finite evidence below separates raw observed target agreement from qualified agreement. Ordinary sector observables use that state's checks, not another sector's residual. Excitations additionally require a qualified lowest AVAILABLE ground reference. Above-gap or degenerate comparison states are excluded individually; partial coverage remains explicit.

| Finite Comparison | Old | New | Scalars | Qualified | Observed Requested | Qualified Requested |
| --- | --- | --- | --- | --- | --- | --- |
| dmrg_vs_full_qp | cosh4-qp-qfull | cosh4-dmrg-best_available | 60 | 60 | 60 | 60 |
| dmrg_vs_full_qp | cosh4-qp-qfull | cosh4-dmrg-highest_qualified | 60 | 60 | 60 | 60 |
| dmrg_vs_qp | cosh10-qp-q6 | cosh10-dmrg-best_available | 60 | 0 | 60 | 0 |
| dmrg_vs_qp | cosh12-qp-q6 | cosh12-dmrg-best_available | 60 | 0 | 60 | 0 |
| dmrg_vs_qp | cosh6-qp-q10 | cosh6-dmrg-best_available | 60 | 15 | 60 | 15 |
| dmrg_vs_qp | cosh8-qp-q8 | cosh8-dmrg-best_available | 60 | 0 | 60 | 0 |
| finite_bath_best | cosh10-qp-q6 | cosh12-qp-q6 | 60 | 60 | 1 | 1 |
| finite_bath_best | cosh4-qp-qfull | cosh6-qp-q10 | 60 | 60 | 1 | 1 |
| finite_bath_best | cosh6-qp-q10 | cosh8-qp-q8 | 60 | 60 | 1 | 1 |
| finite_bath_best | cosh8-qp-q8 | cosh10-qp-q6 | 60 | 60 | 1 | 1 |
| finite_bath_historical | cosh4-qp-qfull | cosh6-qp-q8 | 60 | 60 | 1 | 1 |
| finite_bath_size | cosh2-qp-full | cosh4-qp-qfull | 60 | 45 | 1 | 1 |
| qp_cutoff | cosh10-qp-q4 | cosh10-qp-q6 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh12-qp-q4 | cosh12-qp-q6 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh4-qp-q2 | cosh4-qp-q4 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh4-qp-q4 | cosh4-qp-q6 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh4-qp-q6 | cosh4-qp-qfull | 60 | 60 | 60 | 60 |
| qp_cutoff | cosh6-qp-q4 | cosh6-qp-q6 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh6-qp-q6 | cosh6-qp-q8 | 60 | 60 | 60 | 60 |
| qp_cutoff | cosh6-qp-q8 | cosh6-qp-q10 | 60 | 60 | 60 | 60 |
| qp_cutoff | cosh8-qp-q4 | cosh8-qp-q6 | 60 | 60 | 1 | 1 |
| qp_cutoff | cosh8-qp-q6 | cosh8-qp-q8 | 60 | 60 | 60 | 60 |

| Bath | DMRG Selection | Complete Four States | Qualified | Ground Qualified | Whole Source Calls |
| --- | --- | --- | --- | --- | --- |
| cosh4 | best_available | True | True | True | cosh4-dmrg-even-pair-chi256, cosh4-dmrg-odd_down-chi256, cosh4-dmrg-odd_up-chi256 |
| cosh4 | highest_qualified | True | True | True | cosh4-dmrg-even-pair-chi256, cosh4-dmrg-odd_down-chi256, cosh4-dmrg-odd_up-chi256 |
| cosh6 | best_available | True | False | True | cosh6-dmrg-even-chi1024, cosh6-dmrg-odd_down-chi1024, cosh6-dmrg-odd_up-chi1024 |
| cosh6 | highest_qualified | False | False | False | cosh6-dmrg-odd_down-chi1024 |
| cosh8 | best_available | True | False | False | cosh8-dmrg-even-chi512, cosh8-dmrg-odd_down-chi512, cosh8-dmrg-odd_up-chi512 |
| cosh8 | highest_qualified | False | False | False |  |
| cosh10 | best_available | True | False | False | cosh10-dmrg-even-chi512, cosh10-dmrg-odd_down-chi512, cosh10-dmrg-odd_up-chi512 |
| cosh10 | highest_qualified | False | False | False |  |
| cosh12 | best_available | True | False | False | cosh12-dmrg-even-chi256, cosh12-dmrg-odd_down-chi256, cosh12-dmrg-odd_up-chi256 |
| cosh12 | highest_qualified | False | False | False |  |

### Negative Findings

| Task | Signed Odd Splitting | Exact WS=0 Limit | Qualified |
| --- | --- | --- | --- |
| nrg-z1-k3000 | 0.03168831792 | 0.03 | False |
| nrg-k200 | 0.05285414221 | 0.03 | False |
| nrg-k500 | 0.03748789277 | 0.03 | False |
| nrg-length22 | 0.03748789277 | 0.03 | False |

Historical truncation-invalid points remain invalid even when their disposition allowed later work. Incomplete z grids are not renormalized; complete but invalid grids retain descriptive means marked unqualified. No NRG residual is invented.

Preserved: 157 completed scientific records, 2616 original root rows, plus 2048 Lambda8 ED/QP control roots; 88 compact historical noncompletion outcomes.

| Historical Stage | Outcome | Count |
| --- | --- | --- |
| continuum | deferred | 1 |
| continuum | unstarted | 66 |
| continuum-continue | deferred | 9 |
| continuum-continue | not_needed | 6 |
| continuum-parallel | incomplete | 3 |
| pilot | failed | 1 |
| pilot | unstarted | 2 |

Private-original checksums are reference labels only. Public content-addressed tables and the dataset checksum describe the new payload. Runtime versions and canonical relative source identities are recorded provenance, not executable rehashes. Measured work-duration sums are not parallel wall time or CPU time; RSS/VM are measured peaks.

| Backend | Completed Records | Measured Work Seconds | Peak RSS GiB | Peak VM GiB |
| --- | --- | --- | --- | --- |
| nrg | 69 | 422987.837 | 34.79010773 | 42.19311905 |
| qp | 20 | 33195.85027 | 10.1934433 | 10.77165985 |
| dmrg | 47 | 298633.6255 | 30.63674545 | 30.89536285 |
| ed | 21 | 5.990804657 | 0.01091003418 | 0 |
