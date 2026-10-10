# Standalone Test3 NRG Pilot

**PILOT ONLY, NOT CONVERGED.** Copy this entire directory, including `physics.py`,
outside the repository. It requires Python 3.11+, POSIX shell/flock, and explicit
paths to your NRG Ljubljana tools. Only extraction requires NumPy and h5py.
There are no imports from the installed QD solver or other comparison cases.

The physical default is a single flat reservoir with $\Delta=1$, $D=100$,
$U=2$, $\Gamma=0.4$, detuning $0.5$, and local $B=0.03$. Positive $B$ raises
spin up. The bath is not polarized. Primary states are the lowest even $M=0$
state and odd $M=-1/2,+1/2$ states. No transverse observable is requested.

## Prepare First

Edit `param.loop` **before preparation**. The editable pilot defaults are
`Nz=8`, `Lambda=2`, `Nmax=42`, `keep=2000`, `keepmin=400`, `keepenergy=14`.
These are not production settings or a convergence claim. Start with a single-z
pilot (`#PRELUDE: $Nz = 1;`, giving $z=1$); inspect it before any larger sweep.
The parent study owns all real pilot and subsequent convergence runs.

```sh
export LOOPER=/absolute/path/to/nrgljubljana/scripts/looper
sh ./1_run --prepare-only
```

`1_run` without arguments also only prepares. The controller checks that
`#AUTOLOOP: true` can only generate inputs, then validates every generated deck.
It records `sweep.json` with the template and per-deck SHA-256 hashes. Do not
edit prepared inputs or reuse a failed/partly populated folder. Archive it and
prepare a fresh copy instead. Your tools are never found implicitly through PATH.

## Explicit Execution

Only after approval to run your pilot, set these paths:

```sh
export NRGINIT=/absolute/path/to/nrginit
export NRG=/absolute/path/to/build/c++/nrg
export MATH=/absolute/path/to/math
export NRG_THREADS=1
```

Inside the prepared `1/` directory, `sh ../run_one` explicitly runs the
initializer followed by NRG. In a **fresh, unprepared** copy,
`sh ./1_run --execute` prepares and runs the selected grid sequentially.
The initializer receives the explicit Mathematica path. Logs and per-run scratch
directories are retained. Failed steps stop execution. A completion receipt
hashes `param`, `data`, `raw.h5`, `DONE`, and both logs; it certifies file identity,
not scientific validity. Neither this template's creation nor its offline tests
launch NRG or Mathematica.

## Process

```sh
python3 -m pip install numpy h5py
sh ./2_proc
```

Processing requires the entire prepared sweep and matching completion receipts.
It writes `processed/per_z.csv`, `processed/averages.csv`, and
`processed/summary.json`. The JSON includes all retained states, per-twist
coverage, excluded averages, complete final-sector energies, chain coefficients,
raw energy offsets, and limitations. No reference-mean JSON is supplied: the
study's single-z pilot is not an invented z average.

Only primary IDs present, nondegenerate, and strictly subgap at **every** twist
receive equal-weight means. Extra or missing subgap states are explicitly
flagged. All subgap roots and at least the first two available roots per spin
sector are retained for exclusion diagnostics. Same-sector degeneracy marks
basis-dependent diagonal expectations; it does not abort unrelated states.
Twist min/max are not error bars, and NRG residuals are reported as unavailable.

## Source Conventions

The adapter follows NRG Ljubljana `nrginit/hamiltonian.m` (SIAM/MAGFIELD and
centered bath onsite terms), `nrginit/operators.m` (scalar `SZd`, `SZ2d`),
`nrginit/initial.m` (coefficient units), `c++/invar.hpp` (`SSZ=2M+1`),
`c++/stats.hpp` and `c++/eigen.hpp` (raw energy offsets), and
`c++/truncation.hpp` (zero-based `keepall` indices). Signed SSZ labels include
zero and negative integers. The field belongs in `[extra]`, not `[param]`.
Only `n_d n_d_ud SZd SZ2d` are requested; no `sigma_d` or spectrum is needed.

`render_param(physical, numerical, units='gap')` also accepts `units='bandwidth'`
and optional numerical `untruncated=True` for tiny controls. `read_param(path)`
returns the validated `param` and `extra` sections. `extract_result(directory)`
is read-only and independent of workflow receipts; it returns `physical` in
canonical gap units, `numerical`, `states`, `coverage`, and `diagnostics`.
`diagnostics.deck_gap_scale` records the gap in raw deck units. Energies use

```math
E/\Delta=E_{\mathrm{raw}}/\Delta-E_{\mathrm{bath,vac}}/\Delta-\mathrm{detuning}/\Delta,
\qquad E_{\mathrm{bath,vac}}/\Delta=-\sum_j\sqrt{(\xi_j/\Delta)^2+1}.
```

Here $\xi_j$ are eigenvalues of this twist's normal Wilson chain, not its
hopping coefficients. Its `data` tables supply `Nmax+1` bath sites; the final
unused hopping is removed before diagonalization. Longitudinal bath observables
follow fixed-$M$ identities in local `physics.py`, with physical consistency
checks at absolute tolerance $10^{-8}$.
