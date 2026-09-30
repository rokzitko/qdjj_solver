"""Compatibility for python -m qpsolver.cli."""
from qdjj_solver.qp_solver.cli import main
if __name__ == "__main__":
    raise SystemExit(main())
