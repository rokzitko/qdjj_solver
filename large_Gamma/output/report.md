# Deadline-limited large-hybridization results

**Status:** final deadline-limited report; remaining bath-resolution limits are explicitly reported.

Updated: 2026-09-24 06:15 CEST. Requested completion: 2026-09-24 07:00 +0200.

## Revised accuracy and coverage

The deadline profile targets **$10^{-5}$ finite-bath observable accuracy** and **$10^{-3}$ empirical bath resolution**. The original tighter measurements are retained. Exact ED/BdG references retain their $10^{-8}$ reporting floor.

| Method | Measured base-grid points | Accepted finite-bath observables | All observables bath-resolved |
|---|---:|---:|---:|
| 2QP | 288/288 | 288/288 | 129/288 |
| 4QP | 288/288 | 288/288 | 49/288 |
| 6QP | 288/288 | 288/288 | 55/288 |
| DMRG | 288/288 | 288/288 | 0/288 |

A completed deadline-limited report and a bath-converged calculation are distinct statuses. The last column and `convergence.csv` identify the achieved physical resolution, including cases whose target remains unresolved.

## DMRG accuracy evidence

Where exact same-bath ED/BdG is available, DMRG stops when every reported observable agrees within $10^{-5}$. ED also identifies the lowest eta sector of each parity, so unnecessary excited eta sectors can be omitted from that DMRG run. Degenerate eta minima are retained. This is an observable-specific exact-reference check; the full-state residual and its acceptance flag remain separate diagnostics.

Without an exact reference, both eta signs are solved. Acceptance requires two stable bond refinements, the finest-state sweep check, and a full residual below $10^{-3}\Delta$. Only accepted measurements can define reported breakdown boundaries.

Largest accepted DMRG discrepancy from a same-bath reference in this snapshot: **9.97e-06**, in the normalization of the corresponding observable.

## Signed-gap breakdown brackets

The tables give the **first sampled upward crossing of an absolute signed-gap error of $10^{-3}\Delta$**. Each interaction/phase row uses the largest cosh bath with a complete matched-reference Gamma grid for all three QP cutoffs. These are fixed-bath cutoff-error brackets. Nonmonotonic re-entry and other observables/thresholds are retained in `breakdown.csv`; continuum claims additionally require the stated bath checks.

### Phase difference 0 radians

| U/Delta | Signed levels per lead | 2QP Gamma/Delta | 4QP Gamma/Delta | 6QP Gamma/Delta |
|---:|---:|---|---|---|
| 0 | 16 | 0.1–0.2 | 1–1.5 | 5–6 |
| 0.1 | 10 | 0.1–0.2 | 1–1.5 | 5–6 |
| 0.5 | 10 | 0.05–0.1 | 1.5–2 | 5–6 |
| 1 | 10 | 0.05–0.1 | 1.5–2 | 5–6 |
| 2 | 10 | 0.05–0.1 | 1.5–2 | 5–6 |
| 4 | 10 | 0.02–0.05 | 2–3 | 6–8 |
| 8 | 10 | 0.05–0.1 | 1–1.5 | 6–8 |
| 12 | 10 | 0.05–0.1 | 1–1.5 | 6–8 |

### Phase difference 1.88495559215 radians

| U/Delta | Signed levels per lead | 2QP Gamma/Delta | 4QP Gamma/Delta | 6QP Gamma/Delta |
|---:|---:|---|---|---|
| 0 | 16 | 0.1–0.2 | 1–1.5 | 6–8 |
| 0.1 | 10 | 0.1–0.2 | 1–1.5 | 6–8 |
| 0.5 | 10 | 0.1–0.2 | 1–1.5 | 6–8 |
| 1 | 10 | 0.05–0.1 | 1.5–2 | 6–8 |
| 2 | 10 | 0.02–0.05 | 1.5–2 | 6–8 |
| 4 | 10 | 0.02–0.05 | 1–1.5 | 6–8 |
| 8 | 10 | 0.05–0.1 | 0.75–1 | 6–8 |
| 12 | 10 | 0.05–0.1 | 1–1.5 | 6–8 |

## Limits and retained evidence

Recorded numerical/resource interruptions: **63**. A deadline or resource stop is recorded explicitly and is never treated as convergence.

The complete input and original-run provenance are in `manifest.json`. Scalar checkpoints are retained in the ignored run directories. `observables.csv`, `comparisons.csv`, `convergence.csv`, `bond_convergence.csv` and `failures.json` preserve the measured values and qualification flags. The input files preserve all eight interactions, both phases, and the full initial Gamma grid through 10. No large state files are retained.
