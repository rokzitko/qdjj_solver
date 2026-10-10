# Test4: Model and Parameter Mapping

This document specifies the same physical left/right (L/R) junction in the
ASQ variational manuscript, the canonical NRG `model.m`, and `qdjj_solver`.
The selected representation is **complex, physical L/R**, not the real
even/odd representation of stock NRG SDS models. The default direct bridge is
zero, but the model and current definitions include simultaneous nonzero
$W_N$ and $W_S$ for future comparisons.

This is a model contract, not an execution preset or convergence certificate.
The [final scientific report](README.md) presents the results first; its
[supplement](README.md#supplemental-material) and the
[generated numerical report](output/report.md) retain the matched-model checks
and their limitations. Agreement on small finite systems does not establish
continuum accuracy.

## Sources

- Teodor Ilicin and Rok Zitko, *Quasiparticle-resolved variational theory of
  Andreev spin qubits*, manuscript: the Model section, the definition of the
  Josephson current, and the longitudinal Zeeman term in the magnetic-response
  section. No public DOI or arXiv identifier is recorded in the
  [method references](../../docs/references.md#qp-variational-method).
- The [solver conventions](../../docs/models.md), especially impurity energy
  references, reservoir measure, contact matrices, and physical derivatives.
  The implementation is
  [common/models.py](../../src/qdjj_solver/common/models.py), including
  `reference_model` and `make_model`.
- [NRG Ljubljana](https://github.com/rokzitko/nrgljubljana), upstream modules:
  `nrginit/hamiltonian.m` defines the centered bath onsite terms, pairing
  order, default impurity, and script-model loading;
  `nrginit/initial.m` and `nrginit/wilson.m` define parameter substitution,
  band rescaling, and contact normalization; `nrginit/models.m` contains the
  stock SDS variants; `nrginit/basis.m` defines the SPU1 labels; and
  `c++/h5.hpp` with `c++/operators.hpp` defines complex operator output.
  These are upstream module identifiers, not files shipped by this repository.

Record the actual installed solver and initializer versions for a calculation;
an upstream source checkout need not be identical to the executable in use.
The manuscript's archived NRG settings and tiny spin-selection field are not
test4 defaults or evidence for this parameter point.

## Physical Parameters

All dimensional values in the following table are in units of the positive
physical gap $\Delta=1$; phases are in radians. The bands are identical,
unpolarized, flat, and of finite normal-state half-bandwidth $D=100$, not
full bandwidth 100.
Thus the density of states per spin in each lead is $\rho=1/(2D)=0.005$.

| Quantity | Physical value | Canonical NRG input | `qdjj_solver` input |
|---|---:|---|---|
| Normal half-bandwidth $D$ | 100 | `[param] bandrescale=100` | `bath.bandwidth=100` |
| Positive gap $\Delta$ | 1 | `bcsgap=-0.01` for both leads, before band rescaling | `bath.delta=1` |
| Interaction $U$ | 2 | `[param] U=2` | `u=2` |
| Total hybridization $\Gamma=\Gamma_L+\Gamma_R$ | 0.4 | `[param] Gamma=0.4` | `gamma=0.4` |
| Per-lead hybridization $\Gamma_L=\Gamma_R$ | 0.2 | Obtained by dividing the total hopping amplitude by $\sqrt{2}$ | Split by `geometry="junction"` |
| Preferred dot occupancy $\nu$ | 0.75 | `[extra] nu=0.75` | Mapped to `detuning` |
| Paper detuning $\delta_{\mathrm{paper}}=\nu-1$ | -0.25 | Not the stock NRG `delta` | Dimensionless, not the solver `detuning` |
| Energy detuning $\epsilon=U(1-\nu)$ | 0.5 | Implicit in the canonical dot Hamiltonian | `detuning=0.5` |
| Local Zeeman energy $B_z$ | 0.03 | `[extra] Bz=0.03` | `field=0.03` |
| Phase $\phi=\phi_R-\phi_L$ | $3\pi/5$ | `[extra] phi=1.8849555921538759` | `phi=3*pi/5` |
| Scalar direct hopping $W_N$ | 0 | `[extra] WN=0`, dimensional | `rho_wn=WN/(2*D)` |
| Spin-dependent direct hopping $W_S$ | 0 | `[extra] WS=0`, dimensional | `rho_ws=WS/(2*D)` |

The names have different types: `Gamma` and `gamma` are energies, $V$ and
`gammaPol` below are hopping amplitudes with energy units, $\nu$ and the paper
$\delta$ are dimensionless, and the solver `detuning` is an energy. Bath
`delta` means a positive superconducting gap; stock NRG `delta` means an
impurity energy detuning. None of these uses of `delta` is interchangeable.

## Physical Hamiltonian

Let $d=(d_\uparrow,d_\downarrow)^T$, with spin order `(up, down)`, and define
the normal-measure contact operator in each lead by

```math
C_{l\sigma}=\frac{1}{\sqrt{N}}\sum_k c_{lk\sigma}
\quad\longrightarrow\quad
C_{l\sigma}=\sum_i\sqrt{w_{li}}c_{li\sigma},
\qquad w_{li}\simeq\rho\,d\xi.
```

The finite-bath index $i$ labels a paired orbital, incorporating the paper's
time-reversed partner of $k$. With all Hermitian conjugates included, the
paper-reference Hamiltonian is

```math
\begin{aligned}
H_{\mathrm{paper}}&=H_d^{\mathrm{paper}}+H_L+H_R+H_c+H_T,\\
H_d^{\mathrm{paper}}&=\frac{U}{2}(n_d-\nu)^2+B_z S_d^z,
\qquad S_d^z=\frac{n_\uparrow-n_\downarrow}{2},\\
H_l&=\sum_{i,\sigma}\xi_{li}c_{li\sigma}^\dagger c_{li\sigma}
-\Delta\sum_i\left(c_{li\uparrow}^\dagger c_{li\downarrow}^\dagger
+c_{li\downarrow}c_{li\uparrow}\right),\\
H_c&=\frac{V}{\sqrt{2}}\left[
e^{-i\phi/4}C_L^\dagger d+e^{+i\phi/4}C_R^\dagger d
\right]+\mathrm{h.c.},\\
H_T&=e^{-i\phi/2}C_L^\dagger W_{LR}C_R+\mathrm{h.c.},\\
W_{LR}&=W_N I+iW_S\sigma_z
=\begin{pmatrix}W_N+iW_S&0\\0&W_N-iW_S\end{pmatrix},\\
\Gamma&=\pi\rho V^2,\qquad
V=\sqrt{\frac{2D\Gamma}{\pi}},\qquad
\phi_L=-\phi/2,\quad\phi_R=+\phi/2.
\end{aligned}
```

The coefficient of $C_L^\dagger C_R$ is the displayed matrix, not its
conjugate or its inverse orientation. In `make_model`, tunneling matrices
have **lead-spin rows and dot-spin columns**, and `direct[(0, 1)]` has
L-spin rows and R-spin columns. For this symmetric model
$T_L=T_R=(V/\sqrt{2})I$. The phase factors multiply these matrices as shown;
they are not additional phases of the real pairing amplitudes.

The manuscript also writes
$W=\sqrt{W_N^2+W_S^2}$, $W_N=W\cos(\chi/2)$, and
$W_S=W\sin(\chi/2)$. Thus its spin-resolved bridge coefficient is
$W_\sigma=e^{-i\phi/2}(W_N+i\sigma W_S)$ for $\sigma=+1$ (up) or
$\sigma=-1$ (down). Both bridge inputs here are real, dimensional amplitudes,
not $\rho W$ values. The field is local to the dot, with no lead Zeeman term;
it corresponds to $b_d=B_z$ and $b_{\mathrm{lead}}=0$ in the manuscript.
Positive $B_z$ raises spin up. The auxiliary manuscript generator
$\lambda(n_\uparrow-n_\downarrow)$ would instead have $B_z=2\lambda$.

For discretized contacts, dot hopping includes $\sqrt{w_{li}}$ and bridge
hopping includes $\sqrt{w_{Li}w_{Rj}}$. Do not normalize coarse quadrature or
surrogate weights to sum to one: that changes the physical contacts. The
simple NRG normalization below is restricted to the stated identical flat
bands; it is not a prescription for arbitrary densities, polarized bands,
asymmetric lead couplings, or independently fitted contact normalizations.

## Dot Energy Reference

The canonical `model.m` replaces the default NRG `H1` with the paper dot
Hamiltonian, including its constant. Expanding the square gives

```math
\begin{aligned}
H_d^{\mathrm{paper}}
&=\frac{U}{2}(n_d-1)^2+\epsilon(n_d-1)+B_zS_d^z+C_\nu,\\
\epsilon&=-U\delta_{\mathrm{paper}}=U(1-\nu)=0.5,\\
C_\nu&=\frac{U}{2}(\nu-1)^2=0.0625,\\
H_d^{\mathrm{solver}}&=H_d^{\mathrm{paper}}-C_\nu.
\end{aligned}
```

The conventional Anderson level is $\epsilon_d=\epsilon-U/2=-0.5$.
As an analytic atomic-limit check, the paper energies in the order empty,
up, down, doubly occupied are $0.5625$, $0.0775$, $0.0475$, and $1.5625$.
Subtracting $C_\nu$ gives solver energies $0.5$, $0.015$, $-0.015$, and
$1.5$. These are algebraic reference values, not reported solver results.

Stock NRG instead initializes
$H_{1,\mathrm{stock}}=U(n_d-1)^2/2+\delta_{\mathrm{stock}}n_d$,
before any field addition. Its `delta` is an **energy** and would be
$\delta_{\mathrm{stock}}=\epsilon$, not $\nu-1$.
This stock Hamiltonian differs from the centered solver impurity by
$\epsilon$. That old correction does **not** apply after the canonical
`H1` override: do not subtract 0.5 from canonical test4 energies.

### Two Lead Vacua

Each method must subtract its own **two disconnected BCS lead vacua**,
computed for its actual finite bath or Wilson chain. In the solver's
uncentered electron-bath convention the reference is

```math
E_{\mathrm{BCS},\mathrm{solver}}^{(0)}
=\sum_{l=L,R}\sum_i
\left(\xi_{li}-\sqrt{\xi_{li}^2+\Delta^2}\right).
```

NRG uses onsite terms $\zeta_{ln}(n_{ln}-1)$, not
$\zeta_{ln}n_{ln}$. For the normal-chain eigenvalues $\varepsilon_{la}$
of each retained disconnected chain, with uniform pairing and no anomalous
hopping, its centered reference is

```math
E_{\mathrm{BCS},\mathrm{NRG}}^{(0)}
=-\sum_{l=L,R}\sum_a\sqrt{\varepsilon_{la}^2+\Delta^2}.
```

Use the actual included sites and bonds, not an unused terminal coefficient
or another run's vacuum. Restore the NRG ground-state offset and shell
scaling before treating a saved excitation as an absolute energy. If
$E_{\mathrm{raw}}^{\mathrm{NRG}}$ is that restored energy in the chosen
deck unit, the two useful conventions are

```math
\begin{aligned}
E_{\mathrm{paper}}^{\mathrm{NRG}}
&=E_{\mathrm{raw}}^{\mathrm{NRG}}
-E_{\mathrm{BCS},\mathrm{NRG}}^{(0)},\\
E_{\mathrm{centered}}^{\mathrm{NRG}}
&=E_{\mathrm{paper}}^{\mathrm{NRG}}-C_\nu,\\
E_{\mathrm{paper}}^{\mathrm{solver}}
&=E_{\mathrm{reported}}^{\mathrm{solver}}+C_\nu.
\end{aligned}
```

Compare either the two paper-reference energies or the two centered energies,
never one of each. Divide energies and constants by the physical gap when
reporting gap units. With a nonzero bridge, the reference remains the two
**disconnected** vacua: its physical ground-energy shift must not be removed
by subtracting a bridge-coupled bath vacuum instead. Excitation differences
and the phase derivative are unaffected by $C_\nu$ at fixed $U$ and $\nu$.

## Canonical NRG Mapping

In the seed, `d[]` is the dot, `f[0]` is the L lead's zeroth Wilson site,
and `f[1]` is the R lead's zeroth Wilson site. The latter is a channel label,
not the next site of the L chain. Write these contacts as $f_L$ and $f_R$.
The seed contact terms are exactly

```math
\begin{aligned}
H_{c,\mathrm{seed}}&=\frac{\mathrm{GammaPol}}{\sqrt{2}}
\sum_\sigma\left(
e^{-i\phi/4}f_{L\sigma}^\dagger d_\sigma
+e^{+i\phi/4}f_{R\sigma}^\dagger d_\sigma
\right)+\mathrm{h.c.},\\
H_{T,\mathrm{seed}}&=\sum_\sigma
e^{-i\phi/2}(W_N+i\sigma W_S)
f_{L\sigma}^\dagger f_{R\sigma}+\mathrm{h.c.}
\end{aligned}
```

`GammaPol` in the displayed equation denotes NRG's actual symbol `gammaPol`,
the real total contact amplitude, not another physical input named `Gamma`. For `band=flat` and
`discretization=Z`, the NRG normalization is

```math
\mathrm{GammaPol}=\sqrt{\frac{2\,\mathrm{bandrescale}\,\mathrm{Gamma}}{\pi}}
=V.
```

This follows from `thetaCh=2`,
`theta0Ch=bandrescale*thetaCh*realGamma`, and the unit correction factor
`gammaA=1` for the Z discretization. The implementation uses the documented
`gammaPol` replacement rule, not the legacy `gammaPolch1`/`gammaPolch2`
aliases, which have a different band-rescaling convention.
The division by `Sqrt[2]` is essential: `Gamma=0.4` is total, and using the
undivided amplitude for both leads would double it. Other discretization
conventions require checking their contact correction, not copying this
formula without qualification.

The bridge appears **only between these two zeroth contacts**. Subsequent
Wilson sites form two independent chains with their own normal hoppings and
pairing. Set `rungs=false`; repeating $W_N$ or $W_S$ between sites along the
chains defines a different Hamiltonian.

NRG's pairing coefficient multiplies
$f_\uparrow^\dagger f_\downarrow^\dagger+\mathrm{h.c.}$ with a plus sign.
Therefore the model requires the shared `bcsgap=-0.01` for both physical leads.
Separate `bcsgap1`/`bcsgap2` inputs are rejected by this canonical model, avoiding
ambiguity about the initializer's shared-gap precedence.
After multiplication by `bandrescale=100`, the coefficient is $-1$,
matching the paper's $-\Delta$. The physical gap is the absolute value of
that coefficient, not a negative excitation threshold:

```math
\Delta_{\mathrm{deck}}
 =\mathrm{bandrescale}\,|\mathrm{bcsgap}|>0.
```

### Physical Deck Fragment

This fragment shows physical and representation settings only. It is **not
a runnable or converged execution preset**: retention, chain length, twist
grid, diagonalization, output, and resource settings must be specified
separately. All test4 NRG cases, including atomic, zero-coupling, and
short-chain controls, require $\Lambda\geq4$. This is a test4 restriction,
not a general NRG theorem; do not inherit the smaller values from test3 or
the manuscript's single-channel controls.

```ini
[param]
model=model.m
symtype=SPU1
band=flat
bandrescale=100
discretization=Z
U=2
Gamma=0.4
bcsgap=-0.01
options=COMPLEX SCHUR

[extra]
nu=0.75
Bz=0.03
phi=1.8849555921538759
WN=0
WS=0
```

Use the canonical `model.m` and its `modeloperators.m`, not a stock SDS
variant with this fragment. `U`, `Gamma`, and `bandrescale` belong in
`[param]`; `nu`, `Bz`, `phi`, `WN`, and `WS` belong in `[extra]`. There is
no stock `delta` input to set the canonical dot. Complex arithmetic must
remain enabled even when $W_N=W_S=0$, since finite-phase L/R hopping is
complex.

### Unit Conversion

For a deck using $D$ as its energy unit rather than $\Delta$, divide every
dimensional input by 100: `bandrescale=1`, `U=0.02`, `Gamma=0.004`, and
`Bz=0.0003`, with `WN` and `WS` also divided by 100. Then
$\Delta_{\mathrm{deck}}=0.01$, $\epsilon_{\mathrm{deck}}=0.005$, and
$C_{\nu,\mathrm{deck}}=0.000625$. The dimensionless `nu`, `phi`, and
`bcsgap=-0.01` are unchanged; `bcsgap` inputs are ratios to the
normal half-bandwidth, not already-rescaled gap energies.

Divide output energies and energy-valued current-operator expectations by
0.01 to recover gap units from this $D=1$ deck. The dimensionless inputs
$\rho W_N=W_N/(2D)$ and $\rho W_S=W_S/(2D)$ are invariant under this unit
change. Merely changing `bandrescale` while leaving the other energy inputs
unchanged changes the physical problem.

## Solver Construction

For a declared bath with `delta=1` and `bandwidth=100`, the matching model is:

```python
from math import pi
from qdjj_solver import reference_model

D = 100.0
WN = 0.0
WS = 0.0
h = reference_model(
    bath,
    u=2.0,
    gamma=0.4,
    phi=3*pi/5,
    geometry="junction",
    detuning=0.5,
    field=0.03,
    rho_wn=WN/(2*D),
    rho_ws=WS/(2*D),
    symmetry=False,
    compress=False,
    bath_reference="isolated",
)
```

`bath` must be supplied explicitly; this example chooses no bath quadrature,
QP cutoff, or DMRG bond dimension. Both backends consume this same Hamiltonian.
`symmetry=False` avoids an eta-coordinate reduction, and `compress=False`
retains the full declared bath rather than freezing dark modes. The bath is
still transformed to canonical QP coordinates; physical electron operators
must use `h.physical_operator` when that transformation is needed.
`h.observables["phase_derivative"]` already contains the correct transformed
physical derivative, including the direct bridge.

## Current Operator

Use the unique model-specific operator name `asqJ` in `modeloperators.m`:

```math
J=\frac{\partial H}{\partial\phi}
=\frac{\partial(H_c+H_T)}{\partial\phi},
\qquad
I=\frac{2e}{\hbar}\langle J\rangle,
\qquad e>0.
```

Differentiate the **complete symbolic hopping Hamiltonian before substituting
the numerical parameters**, holding $U$, $\nu$, $B_z$, $\Gamma$, $W_N$, and
$W_S$ fixed. Differentiating after replacing `phi` with a number loses the
observable. In the physical L/R basis its explicit form is

```math
\begin{aligned}
J={}&\frac{V}{\sqrt{2}}\left[
-\frac{i}{4}e^{-i\phi/4}C_L^\dagger d
+\frac{i}{4}e^{+i\phi/4}C_R^\dagger d
\right]+\mathrm{h.c.}\\
&-\frac{i}{2}e^{-i\phi/2}C_L^\dagger W_{LR}C_R
+\mathrm{h.c.}
\end{aligned}
```

Replace the contacts by the zeroth Wilson sites for the NRG seed operator.
The normalization, in any energy unit, is

```math
\frac{I}{2e\Delta/\hbar}
=\frac{\langle\mathrm{asqJ}\rangle}{\Delta_{\mathrm{deck}}}
=\frac{1}{\Delta}\frac{\partial E}{\partial\phi}.
```

The sign is positive with $\phi=\phi_R-\phi_L$, as in the manuscript.
Using units of $e\Delta/\hbar$ instead would double the reported number.
The isolated BCS reference and $C_\nu$ are phase independent and do not
change this derivative. Check expectations against branchwise energy finite
differences; at crossings use the appropriate one-sided branch derivative.

NRG and the solver define this operator in the **same fixed physical basis**,
so there is no additional off-diagonal operator ambiguity from a
phase-dependent gauge transformation. Do not replace it by a derivative of
a phase-dependent diagonalizing basis without the connection terms.
Independent eigenvectors still have arbitrary phases, and degenerate states
require a subspace comparison; compare transition strengths unless their
eigenvector phases have been aligned.

## Spin and Output

`SPU1` conserves total spin projection $M=S_d^z+S_L^z+S_R^z$. Its integer
NRG label is $2M+1$, not a total-spin multiplicity and not directly $2M$:

| State sector | $M$ | NRG label | Solver `twice_sz` | Solver parity |
|---|---:|---:|---:|---:|
| Even, zero projection | 0 | 1 | 0 | 0 |
| Odd, down projection | -1/2 | 0 | -1 | 1 |
| Odd, up projection | +1/2 | 2 | +1 | 1 |

At finite field these labels do not assert exact singlet/doublet SU(2)
multiplets. Both spin-diagonal bridge terms preserve $M$, including when
$W_N$ and $W_S$ are simultaneously nonzero. A transverse dot field or
spin-flip contact does not fit this SPU1 contract. Compare the relevant
sectors before identifying a global ground state; a selected three-branch
inventory alone does not certify subgap completeness.

In complex NRG HDF5 output, an operator block is stored as the real dataset
`matrix` plus its imaginary companion `matrix-imag`. Reconstruct
`block["matrix"][()] + 1j*block["matrix-imag"][()]`; do not discard imaginary
off-diagonal entries.
A Hermitian observable such as `asqJ` has real diagonal expectations, so
check that their imaginary parts vanish within the declared numerical
tolerance. A real diagonal does not imply a real operator matrix.

## Even/Odd Controls

Stock SDS channels are not physical L and R. At $W_N=W_S=0$, identical
leads admit the phase-independent rotation, applied to every paired orbital
or corresponding Wilson site,

```math
c_e=\frac{c_L+c_R}{\sqrt{2}},\qquad
c_o=\frac{i(c_L-c_R)}{\sqrt{2}}.
```

It gives real dot hopping and opposite pairing coefficients:

```math
\begin{aligned}
H_c&=V\cos(\phi/4)C_e^\dagger d
+V\sin(\phi/4)C_o^\dagger d+\mathrm{h.c.},\\
H_{\mathrm{pair}}
&=-\Delta\sum_i c_{ei\uparrow}^\dagger c_{ei\downarrow}^\dagger
+\Delta\sum_i c_{oi\uparrow}^\dagger c_{oi\downarrow}^\dagger
+\mathrm{h.c.}
\end{aligned}
```

Thus the matched even/odd control uses the pairing-coefficient pair
`(-0.01, +0.01)` in the displayed channel order. It must not be copied into
the physical L/R deck, which uses `(-0.01, -0.01)`. The factor $i$ in the
odd channel is what makes the sine hopping real and flips its pairing sign.
Because this rotation is independent of $\phi$, the current transforms by
the same fixed unitary. A nonzero bridge also has to be transformed in full;
the stock bridge-free SDS Hamiltonian is not its replacement.

Before comparing a stock SDS control, match its longitudinal field and
impurity constant as well as its gap signs and hopping. With
$\delta_{\mathrm{stock}}=\epsilon$, its vacuum-subtracted energy obeys
$E_{\mathrm{stock}}=E_{\mathrm{centered}}+\epsilon$; its paper-reference
counterpart is $E_{\mathrm{stock}}-\epsilon+C_\nu$. This is a correction
specific to the **stock** impurity, not to canonical `model.m` output.
Inspect the actual stock variant: source variants that reverse `H0` or
choose another hopping gauge are not interchangeable by name alone.

## Limiting Checks

- At $\Gamma=W_N=W_S=0$, the dot has the analytic atomic energies above,
  the leads are disconnected, and the current is zero.
- At $\Gamma=0$ but nonzero direct hopping, only the dot is disconnected.
  The bridge still joins two superconductors and can carry nonzero Josephson
  current. Do not force `asqJ` to zero or subtract the bridge energy away.
- At $W_N=W_S=0$, the physical L/R and correctly transformed even/odd
  Hamiltonians provide a representation control, not proof of bath convergence.
- Unit-converted decks must preserve energies in gap units, dimensionless
  spin/charge observables, and $I/(2e\Delta/\hbar)$.

The [recorded controls](output/report.md#controls) implement these checks,
including nonzero scalar and spin-dependent bridges. They validate the model
mapping, not production accuracy. Numerical validation must separately address
QP/MPS truncation, NRG retention and twist dependence, and bath convergence.
Every NRG control remains subject to $\Lambda\geq4$.
