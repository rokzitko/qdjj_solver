# Padé chain-expansion reservoirs

`chain_expansion` constructs the ChE reservoirs of D. Bobok, L. Frk, V. Pokorný
and M. Žonda, *Scalable effective models for superconducting nanostructures*,
PRB **112**, 205418 (2025),
[DOI: 10.1103/mxsl-fc96](https://doi.org/10.1103/mxsl-fc96).
It returns a finite reservoir (`DiscreteBath`) in star coordinates, suitable
for both QP diagonalization and DMRG.

```python
from qdjj_solver import chain_expansion, reference_model, solve, Sector

bath = chain_expansion(4, delta=1.0, bandwidth=100.0)
h = reference_model(bath, u=4.0, gamma=2.0, phi=1.0, compress=False)
state = solve(h, cutoff=None, sector=Sector(0, 0), backend="qp")
```

## Parameters and saved bath coefficients

| Parameter | Meaning | Default |
|---|---|---|
| `levels` | Positive integer: signed spinful bath levels, or chain sites, per reservoir | Required |
| `delta` | Positive finite superconducting gap in the Hamiltonian's energy units | `1.0` |
| `bandwidth` | Positive finite normal-state half-bandwidth; see the wide-band convention below | `100.0` |
| `wide_band` | Boolean selecting analytic wide-band ChE coefficients | `False` |

The returned normal levels and weights have exact mirror symmetry. An odd
`levels` count includes a zero normal-energy level. The quantum-impurity
couplings, interaction and phase are supplied to the model, not to the bath
generator. Increasing `levels` changes coefficients throughout the chain.
It is distinct from increasing the QP cutoff or the MPS bond dimension.

The bath's `metadata` dictionary stores `kind="chain-expansion"`, the finite-
or wide-band `target`, `wide_band`, the chain coefficients `chain_h`, the
coefficient quadrature order, the normal-measure error, and the citation.
Reuse the same coefficients with `DiscreteBath.from_record(bath.record())`,
as for other reservoirs.

In a `qdjj-solver` JSON input file, specify the bath as:

```json
{"kind": "chain-expansion", "levels": 4, "delta": 1.0,
 "bandwidth": 100.0, "wide_band": false}
```

## Finite-band construction

Write $x=(\omega/\Delta)^2$ and $\theta_D=\arctan(D/\Delta)$. The bath kernel
$g_D(i\omega)$ defined in the [bath guide](bath_representations.md#the-finite-band-hybridization-target)
has the positive Stieltjes representation

```math
\Delta g_D(i\omega)=\frac{2}{\pi}
\int_0^{\theta_D}\frac{d\theta}{1+x\cos^2\theta}.
```

For $L=2n$, an $n$-point Gaussian quadrature of the positive measure in
$s=\cos^2\theta$ gives the paper's Padé approximant of numerator degree $n-1$
and denominator degree $n$ in $x$. For $L=2n+1$, the corresponding Radau
quadrature fixes one node at $s=1$, giving the zero normal-energy orbital.
Both match the first $L$ expansion coefficients around $x=0$.

Let $q_j$ be the quadrature weight for the measure $`(2/\pi)\,d\theta`$, expressed
in $s=\cos^2\theta$. For an interior node $s_j$, the two normal levels and
their dimensionless residues are

```math
\xi_{j,\pm}=\pm\Delta\sqrt{\frac{1-s_j}{s_j}},\qquad
r_{j,\pm}=\frac{q_j}{2s_j},\qquad
w_{j,\pm}=\pi\rho\Delta r_{j,\pm},\quad \rho=\frac{1}{2D}.
```

The endpoint has $\xi=0$, $r=q$. This produces the same finite-bath
hybridization as the paper's continued fraction without solving an
ill-conditioned Padé coefficient system.

The coefficients are evaluated using positive Gauss–Legendre integration in $\theta$ with
`max(64, 4*levels)` nodes and a two-pass reorthogonalized Jacobi recurrence.
The transformed variable is scaled for narrow bands. Thus coefficient matching
is numerical, to floating-point/integration precision. Tests check analytic
moments, the paper's explicit $L=1,2,3$ coefficients, and direct chain/star
resolvents; the application also compares full many-body chain/star spectra.

## Wide-band construction and normalization

With `wide_band=True`, the coefficients are the paper's analytic Eq. (17):

```math
h_0=L,\qquad h_\ell=\frac{L^2-\ell^2}{4\ell^2-1},\quad 1\leq\ell\lt L.
```

The normal-state chain has zero onsite energies and nearest-neighbor hoppings
$\Delta\sqrt{h_\ell}$. Its coupling to an impurity is
$\sqrt{\Gamma\Delta h_0}$. Diagonalizing the chain gives the star levels and
contact residues. Orthogonality preserves the uniform singlet pairing.

In this mode, **`bandwidth` sets only the density normalization convention**
$\rho=1/(2D)$ used by `DiscreteBath` and the tunneling matrices. It is not a
physical band cutoff in the target function. The weights scale with $\rho$,
so physical combinations such as $w_i/(\pi\rho)$ are independent of this
normalization parameter. `metadata['target']` explicitly says `wide-band`.
Use the default `wide_band=False` for a calculation at a specified finite $D$.
In wide-band mode, the recorded `measure_error=sum(weights)-1` also depends on
the density convention and is not a physical continuum error estimate.

**Do not normalize the weights to sum to one.** ChE matches the low-frequency
kernel. A short chain need not reproduce the normal-state zeroth moment or the
high-frequency tail. Such normalization would change the prescribed model.

## Numerical interpretation

- ChE is a bath approximation; star–chain transformations are coordinate
  changes; QP/MPS truncations approximate the resulting finite many-body problem.
- Bath-size refinement need not be variational or monotonic because the
  Hamiltonian changes with $L$. Compare both odd and even sequences when relevant.
- Ground-state currents, near-edge excitations and small singlet–triplet
  splittings require separate convergence evidence.
- Wide-band ChE and finite-band NRG describe different reservoir limits even
  when their respective many-body calculations are fully converged.

The [worked comparison](../applications/bobok_2025_chain_expansion/README.md)
reproduces selected Fig. 8 current curves and Fig. 17 common-lead spin physics,
including comparisons with GAL, fitted surrogate baths and independent NRG
figure data. The [bath representation guide](bath_representations.md) covers
the other bath families and their measured tradeoffs.
