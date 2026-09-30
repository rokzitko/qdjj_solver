"""Compatibility façade for the earlier qpsolver API."""
import importlib
import sys
from qdjj_solver import __version__ as __version__
from qdjj_solver.qp_solver import *  # noqa: F403
from qdjj_solver.qp_solver import __all__ as __all__

for _name in ("algebra", "baths", "models"):
    sys.modules[f"qpsolver.{_name}"] = importlib.import_module(f"qdjj_solver.common.{_name}")
    globals()[_name] = sys.modules[f"qpsolver.{_name}"]
for _name in ("basis", "solver", "io", "reference", "_core"):
    sys.modules[f"qpsolver.{_name}"] = importlib.import_module(f"qdjj_solver.qp_solver.{_name}")
    globals()[_name] = sys.modules[f"qpsolver.{_name}"]
