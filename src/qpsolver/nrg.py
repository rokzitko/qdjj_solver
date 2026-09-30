"""Compatibility for external-reference comparisons."""
from qdjj_solver.qp_solver.nrg import *  # noqa: F403
from qdjj_solver.qp_solver.nrg import main
if __name__ == "__main__":
    raise SystemExit(main())
