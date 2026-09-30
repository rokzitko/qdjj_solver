"""Resolve and exercise a wheel in a fresh virtualenv away from the checkout."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import venv


CHECK = r'''
import importlib.abc
from pathlib import Path
import sys
target = Path(sys.argv[1]).resolve()
class NoTenpy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] == "tenpy":
            raise ModuleNotFoundError("optional dependency blocked for import test", name="tenpy")
blocker = NoTenpy()
sys.meta_path.insert(0, blocker)
import qdjj_solver
import qpsolver
from qdjj_solver.qp_solver import _core
assert Path(qdjj_solver.__file__).resolve().is_relative_to(target)
assert Path(qpsolver.__file__).resolve().is_relative_to(target)
assert Path(_core.__file__).resolve().is_relative_to(target)
assert qpsolver.Hamiltonian is qdjj_solver.Hamiltonian
h = qdjj_solver.reference_model(qdjj_solver.cosh_grid(1, bandwidth=4), geometry="single")
qp = qdjj_solver.solve(h, options={"eigenpairs": 2})
assert max(qp.residuals) < 1e-9
from qpsolver.io import save_result
save_result(qp, "qp.json", save_vectors=True)
restored = qdjj_solver.load_result("qp.json")
assert (restored.vectors == qp.vectors).all()
qdjj_solver.save_result(qp, "common.json", save_states=True)
assert (qdjj_solver.load_result("common.json").qp_weights == qp.qp_weights).all()
from qpsolver.reference import chain_coefficients
assert len(chain_coefficients(qdjj_solver.cosh_grid(2))[0]) == 4
sys.meta_path.remove(blocker)
if sys.argv[2] == "dmrg":
    from qdjj_solver.dmrg_solver.io import save_result, load_result
    mps = qdjj_solver.solve(h, backend="dmrg", options={"eigenpairs": 2, "require_convergence": True})
    import numpy as np
    np.testing.assert_allclose(mps.energies, qp.energies, atol=2e-10, rtol=0)
    save_result(mps, "mps.json", save_states=True)
    restored = load_result("mps.json")
    np.testing.assert_allclose(abs(restored.overlaps(mps)), np.eye(2), atol=1e-8, rtol=0)
from qdjj_solver.cli import main
try:
    main(["--version"])
except SystemExit as error:
    assert error.code == 0
else:
    raise AssertionError("version command did not exit")
'''


def check_installed(python, target, directory, environment, dmrg):
    subprocess.run([str(python), "-I", "-m", "pip", "check"],
                   cwd=directory, env=environment, check=True)
    subprocess.run([str(python), "-I", "-B", "-c", CHECK, str(target), "dmrg" if dmrg else "qp"],
                   cwd=directory, env=environment, check=True)
    for name in ("qdjj-solver", "qpsolver"):
        executable = python.parent/(name + (".exe" if os.name == "nt" else ""))
        subprocess.run([str(executable), "--version"], cwd=directory, env=environment, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path, nargs="?")
    parser.add_argument("--installed", action="store_true", help="check a wheel already installed by cibuildwheel")
    parser.add_argument("--dmrg", action="store_true", help="install and check the DMRG/HDF5 extra")
    args = parser.parse_args()
    if (args.wheel is None) == (not args.installed):
        parser.error("provide a wheel or --installed")
    wheel = None if args.wheel is None else args.wheel.resolve()
    if wheel is not None and (not wheel.is_file() or wheel.suffix != ".whl"):
        parser.error("provide one built wheel")
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    for variable in ("PYTHONPATH", "PYTHONHOME"):
        environment.pop(variable, None)
    with tempfile.TemporaryDirectory(prefix="qdjj-wheel-") as directory:
        if args.installed:
            check_installed(Path(sys.executable), Path(sys.prefix), directory, environment, args.dmrg)
            return
        target = Path(directory)/"venv"
        venv.EnvBuilder(with_pip=True).create(target)
        scripts = target/("Scripts" if os.name == "nt" else "bin")
        python = scripts/("python.exe" if os.name == "nt" else "python")
        requirement = str(wheel) + ("[dmrg]" if args.dmrg else "")
        subprocess.run([str(python), "-I", "-m", "pip", "install", requirement],
                       cwd=directory, env=environment, check=True)
        check_installed(python, target, directory, environment, args.dmrg)
    print("Clean wheel dependencies, imports, native calculation, console scripts, and serialization passed.")


if __name__ == "__main__":
    main()
