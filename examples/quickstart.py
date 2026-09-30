"""Small QP-only calculation with a reproducible finite-model check."""

from qdjj_solver import Sector, cosh_grid, reference_model, solve


bath = cosh_grid(pairs=1, delta=1.0, bandwidth=5.0)
h = reference_model(bath, u=1.6, gamma=0.3, phi=0.8, compress=False)
sector = Sector(parity=1, twice_sz=1)
qp = solve(h, cutoff=2, sector=sector, backend="qp")
print(qp.energies, qp.qp_weights)

assert abs(qp.energies[0] - (-0.2950595967187083)) < 1e-10
assert max(qp.residuals) < 1e-9
assert abs(qp.qp_weights[0].sum() - 1) < 1e-12
